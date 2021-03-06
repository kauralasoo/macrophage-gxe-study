library("devtools")
library("qvalue")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import permutation p-values
naive_qtls = importFastQTLTable("results/SL1344/fastqtl/output/naive_500kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_name_map)
IFNg_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_500kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_name_map)
SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/SL1344_500kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_name_map)
IFNg_SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_SL1344_500kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_name_map)

#Combine all hits to a list
qtl_list = list(naive = naive_qtls, IFNg = IFNg_qtls, SL1344 = SL1344_qtls, IFNg_SL1344 = IFNg_SL1344_qtls)
fastqtl_hits = purrr::map(qtl_list, ~dplyr::filter(.,qvalue < 0.1))
saveRDS(qtl_list, "results/SL1344/eQTLs/fastqtl_min_pvalues.rds")

#Calculate Pi1 statistic
pi1_matrix = calculatePairwisePi1(qtl_list)
write.table(pi1_matrix, "results/SL1344/eQTLs/properties/pi1_matrix.txt", quote = FALSE, sep = "\t")
pi1_matrix_tidy = calculatePairwisePi1(qtl_list, tidy = TRUE)
write.table(pi1_matrix_tidy, "results/SL1344/eQTLs/properties/pi1_matrix_tidy.txt", quote = FALSE, sep = "\t")

#Make a QQ plot
ggplot(naive_qtls, aes(x = p_expected, y = p_beta_log10)) +
  geom_point() +
  stat_abline(slope = 1, intercept = 0, color = "red")


#Import permutation p-values from 100kb window
naive_qtls = importFastQTLTable("results/SL1344/fastqtl/output/naive_100kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_name_map)
IFNg_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_100kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_name_map)
SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/SL1344_100kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_name_map)
IFNg_SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_SL1344_100kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_name_map)
qtl_list = list(naive = naive_qtls, IFNg = IFNg_qtls, SL1344 = SL1344_qtls, IFNg_SL1344 = IFNg_SL1344_qtls)
saveRDS(qtl_list, "results/SL1344/eQTLs/fastqtl_min_pvalues_100kb.rds")

#Import min p-values from fastqtl 100kb and 500kb windows centred around the TSS
fastqtl_100kb_list = list(
  naive = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_centre/naive_100kb_permuted.txt.gz",
  IFNg = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_centre/IFNg_100kb_permuted.txt.gz",
  SL1344 = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_centre/SL1344_100kb_permuted.txt.gz",
  IFNg_SL1344 = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_centre/IFNg_SL1344_100kb_permuted.txt.gz")
fastqtl_100kb_min_pvalues = purrr::map(fastqtl_100kb_list, ~importFastQTLTable(.) %>% enrichFastQTLPvalues(gene_name_map))
saveRDS(fastqtl_100kb_min_pvalues, "results/SL1344/eQTLs/fastqtl_min_pvalues_100kb.rds")

fastqtl_500kb_list = list(
  naive = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_centre/naive_500kb_permuted.txt.gz",
  IFNg = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_centre/IFNg_500kb_permuted.txt.gz",
  SL1344 = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_centre/SL1344_500kb_permuted.txt.gz",
  IFNg_SL1344 = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_centre/IFNg_SL1344_500kb_permuted.txt.gz")
fastqtl_500kb_min_pvalues = purrr::map(fastqtl_500kb_list, ~importFastQTLTable(.) %>% enrichFastQTLPvalues(gene_name_map))
saveRDS(fastqtl_500kb_min_pvalues, "results/SL1344/eQTLs/fastqtl_min_pvalues_500kb.rds")

fastqtl_100kb_n42_list = list(
  naive = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_42/naive_100kb_permuted.txt.gz",
  IFNg = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_42/IFNg_100kb_permuted.txt.gz",
  SL1344 = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_42/SL1344_100kb_permuted.txt.gz",
  IFNg_SL1344 = "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_42/IFNg_SL1344_100kb_permuted.txt.gz")
fastqtl_100kb_42_min_pvalues = purrr::map(fastqtl_100kb_n42_list, ~importFastQTLTable(.) %>% enrichFastQTLPvalues(gene_name_map))
saveRDS(fastqtl_100kb_42_min_pvalues, "results/SL1344/eQTLs/fastqtl_min_pvalues_100kb_42.rds")


