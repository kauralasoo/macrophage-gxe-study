library("ggplot2")
library("gplots")
library("devtools")
load_all("../seqUtils/")

#Load expresison dataset preapred previously by processExpressionData.R script
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata

#Filter the design matrix
design = dplyr::filter(expression_dataset$design, !(donor == "fpdj")) %>% tbl_df() %>% #Remove all fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) %>% #Remove second ougl sample (dium)
  dplyr::filter(!(donor == "mijn")) %>% #Remove mijn (wrong line from CGAP)
  dplyr::filter(!(donor == "jorr")) #Very strong outlier in PEER analysis
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#Keep only one condition
line_meta = dplyr::filter(sample_meta, condition == "A")

#Import PEER results
peer_factors = read.table("results/SL1344/PEER/25_factors/factors.txt", sep =",")
peer_factors = peer_factors[,2:ncol(peer_factors)]
colnames(peer_factors) = paste("PEER_factor_", c(1:ncol(peer_factors)), sep = "")

#Add PEER results to the sample metadata
peer_factors_2 = dplyr::mutate(peer_factors, sample_id = line_meta$sample_id)
expanded_meta = dplyr::left_join(line_meta, peer_factors_2, by = "sample_id")
expanded_meta = dplyr::mutate(expanded_meta, RNAseq_auto = ifelse(chemistry == "V4_auto", "yes", "no"))

#Import precision
precisions = read.table("results/SL1344/PEER/25_factors/precision.txt", sep = ",")
plot(1/precisions$V1)

#Explore the correlation between PEER factors
heatmap.2(cor(peer_factors), margins = c(10,10))

### Identify which factors are most correlated with PEER factors ####
covariates = dplyr::select(expanded_meta, mean_purity_filtered, ng_ul_mean, RNAseq_auto, diff_days, passage_diff)

peer_explained = explainPEER(peer_factors[,1:5], covariates)
peer_explained$peer_factor = factor(peer_explained$peer_factor, 
                                    levels = paste("PEER_factor_", c(1:ncol(peer_factors)), sep = ""))
ggplot(peer_explained, aes(x = peer_factor, y = covariate, fill = sqrt(r_squared))) + 
  geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = scales::muted("blue"), 
                       high = scales::muted("red"), midpoint = 0.4, name = "Correlation")

#Make the same plot by excluding some of the purity outliers
filter = expanded_meta$mean_purity_filtered > 0.95
peer_explained = explainPEER(peer_factors[filter,1:5], covariates[filter,])
peer_explained$peer_factor = factor(peer_explained$peer_factor, 
                                    levels = paste("PEER_factor_", c(1:ncol(peer_factors)), sep = ""))
ggplot(peer_explained, aes(x = peer_factor, y = covariate, fill = sqrt(r_squared))) + 
  geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = scales::muted("blue"), 
                       high = scales::muted("red"), midpoint = 0.4, name = "Correlation")

#Visualize the RNA-Seq batch effect
ggplot(expanded_meta, aes(x = PEER_factor_1, y = PEER_factor_2, label = sample_id, color = RNAseq_auto)) + 
  geom_point() + geom_text()

#Visualize the purity effect with one outlier
ggplot(expanded_meta, aes(x = mean_purity_filtered, y = PEER_factor_1, label = sample_id)) + 
  geom_point() + geom_text()

#Remove bubh outlier:
ggplot(dplyr::filter(expanded_meta, mean_purity_filtered > 0.9), aes(x = mean_purity_filtered, y = PEER_factor_1, label = sample_id)) + 
  geom_point() + geom_text()

#Correlation between factor 1 and RNA concentration
ggplot(expanded_meta, aes(x = ng_ul_mean, y = PEER_factor_1, label = sample_id)) + geom_point() + geom_text()
cor.test(as.numeric(factor(expanded_meta$ng_ul_mean)), expanded_meta$PEER_factor_1, method = "pearson")

