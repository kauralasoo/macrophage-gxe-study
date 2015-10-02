---
title: "Macrophage GxE study QC report"
output: html_document
---


This is a QC report for the macrophage GxE study.

Load required packages and import raw metadata from disk.

```r
setwd("../../")
library("dplyr")
library("ggplot2")
library("tidyr")

dat = tbl_df(read.csv("macrophage-gxe-study/data/sample_lists/line_metadata_020315.csv", stringsAsFactors = FALSE, na.strings = ""))
flow_purity = readRDS("results/covariates/flow_cytometry_purity.rds")
rna_concentrations = readRDS("results/covariates/rna_concentrations.rds")
```

##Load HipSci genotype sample names

```r
genotypes = read.table("/Users/alasoo/projektid/macrophage-gxe-study/genotypes/genotype_sample_names.txt", stringsAsFactors = FALSE)
colnames(genotypes) = c("genotype_id")
genotypes = tidyr::separate(genotypes, genotype_id, into = c("batch_id", "line_id"), sep = "-", remove = FALSE) %>% 
  dplyr::select(genotype_id, line_id)
```

Do some munging of the data before we can start making pretty plots.

```r
#Convert all date fields to date
line_metadata = dat %>%
  tidyr::separate(line_id, into = c("donor","clone"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(ips_received = as.Date(ips_received, "%d/%m/%Y"), 
                ips_started = as.Date(ips_started, "%d/%m/%Y"),
                EB_formation = as.Date(EB_formation, "%d/%m/%Y"),
                diff_start = as.Date(diff_start, "%d/%m/%Y"),
                MF_harvest = as.Date(MF_harvest, "%d/%m/%Y"),
                salmonella = as.Date(salmonella, "%d/%m/%Y"),
                flow_date = as.Date(flow_date, "%d/%m/%Y"),
                terminated = as.Date(terminated, "%d/%m/%Y"),
                extraction_date = as.Date(extraction_date, "%d/%m/%Y"),
                rna_submit = as.Date(rna_submit, "%d/%m/%Y"))

#Add mean RNA concentration and flow purity to the data
mean_rna_concentrations = dplyr::select(rna_concentrations, donor, ng_ul_mean, replicate) %>% unique()
mean_flow_purity = flow_purity %>%
  dplyr::filter(!(donor == "gedo" & channel == "Pacific.Blue.A")) %>% #Remove an outlier measurement
  group_by(donor,flow_date) %>% 
  dplyr::summarize(mean_purity = mean(purity), max_purity = max(purity))

line_data = dplyr::left_join(line_metadata, mean_rna_concentrations, by = c("donor", "replicate")) %>%
  dplyr::left_join(mean_flow_purity, by = c("donor", "flow_date")) %>%
  dplyr::left_join(genotypes, by = "line_id")
```

## Export genotypes for RNA-Seq samples

```r
rna_seq_genotypes = dplyr::filter(line_data, !is.na(rna_submit)) %>% 
  dplyr::select(line_id, genotype_id) %>% 
  arrange(line_id) %>% 
  unique()
write.table(rna_seq_genotypes, "genotypes/rna_seq_genotypes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

```
## Warning in file(file, ifelse(append, "a", "w")): cannot open file
## 'genotypes/rna_seq_genotypes.txt': No such file or directory
```

```
## Error in file(file, ifelse(append, "a", "w")): cannot open the connection
```


#Basic statistics
##Duration
How long does it take from EB formation to stimulation experiment?

```r
duration_df = line_data %>%
  dplyr::filter(status == "Success") %>%
  dplyr::mutate(diff_duration = as.numeric(salmonella - EB_formation)) %>%
  dplyr::mutate(ips_culture_duration = as.numeric(EB_formation - ips_started))
ggplot(duration_df, aes(x = diff_duration)) + 
  geom_histogram(binwidth = 1) +
  xlab("Duration of differentiation (days)")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

How long does it take to expand iPS cells before differentiation?

```r
ggplot(duration_df, aes(x = ips_culture_duration, fill = received_as)) + 
  geom_histogram(binwidth = 1) + 
  xlab("Duration of iPS cell culture (days).")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 

How many medium changes with cytokines (IL-3 + M-CSF) do we need to perform for each line?

```r
ggplot(duration_df, aes(x = medium_changes)) + 
  geom_histogram(binwidth = 1) + 
  xlab("Median number of medium changes per line.")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

From this plot, the median number of medium changes per line is 8. TODO: add cost estimates based on this.

##Success rate
What is the overall success rate of our macrophage differentiations?

```r
success_df = line_data %>%
  dplyr::mutate(success = ifelse(status == "Success", "Success", "Fail")) %>%
  dplyr::mutate(status = ifelse(status == "RNA_QC_fail", "Stimulation_fail", status))
ggplot(success_df, aes(x = success, fill = status)) + geom_bar()
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

Our success rate thus far is 39/59 = 66%. Most commonly, the iPS line fails to differentiate completely and we do not get any macrophages (Diff fail). However, often the line does differentiate but the resulting macrophages are contaminate by some other cell type that cannot be easily separate (FC_QC_fail).

How are these successes and failures distributed over time?


```r
monthly_df = success_df %>%
  dplyr::select(donor, clone, EB_formation, status) %>% 
  dplyr::arrange(EB_formation) %>% 
  dplyr::mutate(month = format(EB_formation, "%b,%y")) %>%
  dplyr::mutate(month = factor(month, levels = unique(month)))
ggplot(monthly_df, aes(x = month, fill = status)) + geom_bar() + 
  xlab("Start of differentiation (month)")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 

We spent first couple of months of 2014 setting everything up and optimizing the stimulation assays. From May-July we were limited by the number of lines that CGaP could provide us. We received our first larger batch of cells in August, but unfortunately many of these differentiations failed, probably because of too infrequent medium change. After that we increased the medium change frequncy to 4-5 days and this improved our success rates dramatically. We did not start any new differentiations after the beginning of November, because we wanted to finish all of them before Christmas.

##Batch sizes
How large are the stimulation batches?

```r
batch_df = line_data %>%
  dplyr::filter(status == "Success") %>%
  group_by(salmonella) %>% 
  summarize(batch_size = length(salmonella))
ggplot(batch_df, aes(x = batch_size)) + 
  geom_histogram(binwidth = 1) +
  xlab("Size of stimulation batch (lines)")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 

The median is 3 lines per stimulation batch, but there are also two batches that only contained 1 line.

# Flow cytometry
## Quantifying the purity of macrophage samples
For all of the iPSDM samples we performed flow cytometry for 3 cell surface markers: CD14, CD16 and CD206. We then counted the percentage of the cells that stained positive for each of marker. Finally, we used the maximum of the three percentage as the purity score of the sample (we found maximum to be more robust than the mean).

Extract data and remove sample in which flow cytometry was done more than 20 days later than RNA extraction.

```r
flow_df = line_data %>% 
  dplyr::filter(status %in% c("Success", "FC_QC_fail"), !is.na(max_purity)) %>%
  dplyr::filter(flow_date - salmonella < 20 | is.na(salmonella)) 
```
Make a histogram of the purity scores. Samples that failed flow cytometry QC are shown in red.

```r
ggplot(flow_df, aes(x = max_purity-0.001, fill = status)) + 
  geom_histogram(binwidth = 0.01) + 
  xlab("Maximum purity")
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png) 

## Outlier samples on PCA
We found that some of the samples on the gene expression PCA plot clearly clustered separately from the others. Do these samples also have lower purity scores? Yes, it does seem to be the case:

```r
cluster_df = flow_df %>%
  dplyr::mutate(separate_cluster = ifelse(donor %in% c("iasn","ougl","debk","gomv","ffdp","peop","huls"), "Yes", "No")) %>%
  dplyr::filter(!is.na(rna_submit)) #Sequenced samples only
ggplot(cluster_df, aes(x = separate_cluster, y = max_purity, label = donor)) +
  geom_boxplot() +
  geom_text() +
  xlab("Is in separate cluster?") + 
  ylab("Maximum purity")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png) 

Furthermore this difference is not explained by RNA concentration:

```r
ggplot(cluster_df, aes(x = separate_cluster, y = ng_ul_mean)) +
  geom_violin() +
  xlab("Is in separate cluster?") + 
  ylab("RNA concentration (ng/ul)")
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png) 