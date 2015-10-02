# Macrophage GxE study QC report
  

```
## Warning: package 'knitr' was built under R version 3.2.2
```

```
## Loading seqUtils
## Creating a generic function for 'nchar' from package 'base' in package 'S4Vectors'
```

Load expresison dataset preapred previously by processExpressionData.R script,

```r
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
```

Filter the design matrix,

```r
design = dplyr::filter(expression_dataset$design, !(donor == "fpdj")) %>% tbl_df() %>% #Remove all fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) %>% #Remove second ougl sample (dium)
  dplyr::filter(!(donor == "mijn")) %>% #Remove mijn (wrong line from CGAP)
  dplyr::filter(!(donor == "jorr")) #Very strong outlier in PEER analysis
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))
line_meta = dplyr::filter(sample_meta, condition == "A")
```

Import PEER results

```r
peer_factors = read.table("results/SL1344/PEER/25_factors/factors.txt", sep =",")
peer_factors = peer_factors[,2:ncol(peer_factors)]
colnames(peer_factors) = paste("PEER_factor_", c(1:ncol(peer_factors)), sep = "")
```

Add PEER results to the sample metadata

```r
peer_factors_2 = dplyr::mutate(peer_factors, sample_id = line_meta$sample_id)
expanded_meta = dplyr::left_join(line_meta, peer_factors_2, by = "sample_id")
expanded_meta = dplyr::mutate(expanded_meta, RNAseq_auto = ifelse(chemistry == "V4_auto", "yes", "no"))
```

Look at the variation explained by each PEER factor,

```r
precisions = read.table("results/SL1344/PEER/25_factors/precision.txt", sep = ",")
plot(1/precisions$V1)
```

![](explorePeerFactors_files/figure-html/unnamed-chunk-6-1.png) 

Explore the correlation between PEER factors

```r
heatmap.2(cor(peer_factors), margins = c(10,10))
```

![](explorePeerFactors_files/figure-html/unnamed-chunk-7-1.png) 

Identify which factors are most correlated with PEER factors

```r
covariates = dplyr::select(expanded_meta, mean_purity_filtered, ng_ul_mean, RNAseq_auto, diff_days, passage_diff)

peer_explained = explainPEER(peer_factors[,1:5], covariates)
peer_explained$peer_factor = factor(peer_explained$peer_factor, 
                                    levels = paste("PEER_factor_", c(1:ncol(peer_factors)), sep = ""))
ggplot(peer_explained, aes(x = peer_factor, y = covariate, fill = sqrt(r_squared))) + 
  geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = scales::muted("blue"), 
                       high = scales::muted("red"), midpoint = 0.4, name = "Correlation")
```

![](explorePeerFactors_files/figure-html/unnamed-chunk-8-1.png) 

Make the same plot by excluding some of the purity outliers

```r
filter = expanded_meta$mean_purity_filtered > 0.95
peer_explained = explainPEER(peer_factors[filter,1:5], covariates[filter,])
peer_explained$peer_factor = factor(peer_explained$peer_factor, 
                                    levels = paste("PEER_factor_", c(1:ncol(peer_factors)), sep = ""))
ggplot(peer_explained, aes(x = peer_factor, y = covariate, fill = sqrt(r_squared))) + 
  geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = scales::muted("blue"), 
                       high = scales::muted("red"), midpoint = 0.4, name = "Correlation")
```

![](explorePeerFactors_files/figure-html/unnamed-chunk-9-1.png) 

Visualize the RNA-Seq batch effect,

```r
ggplot(expanded_meta, aes(x = PEER_factor_1, y = PEER_factor_2, label = sample_id, color = RNAseq_auto)) + 
  geom_point() + geom_text()
```

![](explorePeerFactors_files/figure-html/unnamed-chunk-10-1.png) 

Visualize the purity effect with one outlier,

```r
ggplot(expanded_meta, aes(x = mean_purity_filtered, y = PEER_factor_1, label = sample_id)) + 
  geom_point() + geom_text()
```

```
## Warning: Removed 5 rows containing missing values (geom_point).
```

```
## Warning: Removed 5 rows containing missing values (geom_text).
```

![](explorePeerFactors_files/figure-html/unnamed-chunk-11-1.png) 

Remove bubh outlier,

```r
ggplot(dplyr::filter(expanded_meta, mean_purity_filtered > 0.9), aes(x = mean_purity_filtered, y = PEER_factor_1, label = sample_id)) + 
  geom_point() + geom_text()
```

![](explorePeerFactors_files/figure-html/unnamed-chunk-12-1.png) 

Correlation between factor 1 and RNA concentration

```r
ggplot(expanded_meta, aes(x = ng_ul_mean, y = PEER_factor_1, label = sample_id)) + geom_point() + geom_text()
```

![](explorePeerFactors_files/figure-html/unnamed-chunk-13-1.png) 

```r
cor.test(as.numeric(factor(expanded_meta$ng_ul_mean)), expanded_meta$PEER_factor_1, method = "pearson")
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  as.numeric(factor(expanded_meta$ng_ul_mean)) and expanded_meta$PEER_factor_1
## t = -2.237, df = 55, p-value = 0.02936
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.51090926 -0.03051533
## sample estimates:
##        cor 
## -0.2887868
```