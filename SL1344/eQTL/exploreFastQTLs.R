library("devtools")
library("qvalue")
library("dplyr")
library("ggplot2")
library("SNPRelate")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import genotypes
#SNPRelate::snpgdsVCF2GDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz", 
#                         "results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.gds", method = "copy.num.of.ref")
#vcf_file = gdsToMatrix("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.gds")
#saveRDS(vcf_file, "results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")
#vcf_file = readRDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")

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

