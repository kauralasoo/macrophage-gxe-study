library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("ggplot2")
library("purrr")

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
gene_name_map = dplyr::select(atac_list$gene_metadata, gene_id, gene_name)

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#import IRF2 peaks
irf2_peaks = readRDS("results/public_chromatin/joint_peaks/IRF2_Fairfax_peaks.rds")

#Find ATAC peaks that overlap with IRF2 peaks
atac_peaks = dplyr::transmute(atac_list$gene_metadata, seqnames = chr, start, end, strand, gene_id) %>% 
  dataFrameToGRanges()
irf2_atac_peaks = atac_peaks[GenomicRanges::queryHits(GenomicRanges::findOverlaps(atac_peaks, irf2_peaks)),]$gene_id

#Calculate mean accessibility at IRF2 peaks
irf2_means = colMeans(atac_list$cqn[irf2_atac_peaks,])
mean_accessibility = tidyVector(irf2_means, value_id = "mean_access")

#Etract genotype data
gt_df = vcf_file$genotypes["rs13149699",] %>% 
  tidyVector(sample_id = "genotype_id", value_id = "rs13149699")

#Merge all together
df = dplyr::left_join(atac_list$sample_metadata, gt_df, by = "genotype_id") %>% 
  dplyr::left_join(mean_accessibility, by = "sample_id")
ggplot(df, aes(x = as.factor(rs13149699), y = mean_access)) + 
  facet_wrap(~condition_name) + geom_boxplot() + geom_point()

#Test for association
df_filtered = dplyr::filter(df, condition_name == "IFNg_SL1344")
df_filtered$assigned_millions = df_filtered$assigned/1e6
m1 = lm(mean_access ~ sex_binary + assigned_frac, df_filtered)
m2 = lm(mean_access ~ sex_binary + assigned_frac + rs13149699, df_filtered)
anova(m2,m1)


