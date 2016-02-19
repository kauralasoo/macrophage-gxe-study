library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("Rsamtools")

#Import data
acldl_list = readRDS("results/acLDL/acLDL_combined_expression_data_covariates.rds")
gene_name_map = acldl_list$gene_metadata %>% dplyr::select(gene_id, gene_name, gene_biotype)

#Import SNP coordinates
snp_coords = readr::read_delim("genotypes/acLDL/imputed_20151005/imputed.70_samples.snp_coords.txt", 
                               delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))

#Extract minimal p-value for each condition
ctrl_eigen_pvalue = eigenMTImportResults("results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.eigenMT.txt")
acldl_eigen_pvalue = eigenMTImportResults("results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.eigenMT.txt")

min_pvalue_list = list(Ctrl = ctrl_eigen_pvalue,
                       AcLDL = acldl_eigen_pvalue)
saveRDS(min_pvalue_list, "results/acLDL/eQTLs/acLDL_rasqual_min_pvalues.rds")


#Import fastQTL p-values from disk
ctrl_perm_pvalue = importFastQTLTable("results/acLDL/fastqtl/output/Ctrl_permuted.txt.gz")
acldl_perm_pvalue = importFastQTLTable("results/acLDL/fastqtl/output/AcLDL_permuted.txt.gz")
perm_pvalue_list = list(Ctrl = ctrl_perm_pvalue, AcLDL = acldl_perm_pvalue)
saveRDS(perm_pvalue_list, "results/acLDL/eQTLs/acLDL_fastqtl_min_pvalues.rds")

#Compare the correlation betweeb rasqual and fastQTL p-values
fastqtl_vs_rasqual = ggplot(ctrl_matched, aes(x = -log(p_linear, 10), y = -log(p_rasqual,10))) + geom_point()
ggsave("results/acLDL/eQTLs/fastqtl_vs_rasqual.pdf", plot = fastqtl_vs_rasqual)

#Export lead-qQTLs
dplyr::left_join(ctrl_eigen_pvalue, gene_name_map, by = "gene_id") %>%
  write.table("results/acLDL/eQTLs/ctrl_rasqual_hits.txt", quote = FALSE, row.names = FALSE, sep = "\t")
dplyr::left_join(acldl_eigen_pvalue, gene_name_map, by = "gene_id") %>%
  write.table("results/acLDL/eQTLs/acLDL_rasqual_hits.txt", quote = FALSE, row.name = FALSE, sep = "\t")


