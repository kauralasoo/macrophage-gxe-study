library("dplyr")
library("devtools")
library("Rsamtools")
library("purrr")
load_all("../seqUtils/")

##### Import RASQUAL results #####
#Nominal RASQUAL p-values
ctrl_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/acLDL/rasqual/Ctrl_500kb.eigenMT.txt")
acLDL_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/acLDL/rasqual/AcLDL_500kb.eigenMT.txt")
min_pvalue_list = list(Ctrl = ctrl_eigenMT, AcLDL = acLDL_eigenMT)

#Permutation Rasqual p-values
ctrl_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/acLDL/rasqual/random_permutation/Ctrl_500kb.eigenMT.txt")
acLDL_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/acLDL/rasqual/random_permutation/AcLDL_500kb.eigenMT.txt")
min_pvalue_list_random = list(Ctrl = ctrl_eigenMT, AcLDL = acLDL_eigenMT)

#Calculate empirical FDR threshold
rasqual_fdr_thres_list = purrr::map2(min_pvalue_list, min_pvalue_list_random, 
                                     ~dplyr::mutate(.x, fdr_thresh = getFDR(p_eigen, .y$p_eigen, alpha = 0.1)))
saveRDS(rasqual_fdr_thres_list, "results/acLDL/eQTLs/rasqual_min_pvalues.rds")


#Import data
acldl_list = readRDS("results/acLDL/acLDL_combined_expression_data_covariates.rds")
gene_name_map = acldl_list$gene_metadata %>% dplyr::select(gene_id, gene_name, gene_biotype)

#Import SNP coordinates
snp_coords = readr::read_delim("genotypes/acLDL/imputed_20151005/imputed.70_samples.snp_coords.txt", 
                               delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))

#Extract minimal p-value for each condition
ctrl_eigen_pvalue = eigenMTImportResults("results/acLDL/rasqual/Ctrl_500kb.eigenMT.txt")
acldl_eigen_pvalue = eigenMTImportResults("results/acLDL/rasqual/AcLDL_500kb.eigenMT.txt")

min_pvalue_list = list(Ctrl = ctrl_eigen_pvalue,
                       AcLDL = acldl_eigen_pvalue)
saveRDS(min_pvalue_list, "results/acLDL/eQTLs/acLDL_rasqual_min_pvalues.rds")

#Extract minimal p-value for each condition (filtered)
ctrl_eigen_pvalue = eigenMTImportResults("results/acLDL/rasqual/output_filtered/Ctrl_500kb.eigenMT.txt")
acldl_eigen_pvalue = eigenMTImportResults("results/acLDL/rasqual/output_filtered/AcLDL_500kb.eigenMT.txt")

min_pvalue_list_filtered = list(Ctrl = ctrl_eigen_pvalue,
                       AcLDL = acldl_eigen_pvalue)
saveRDS(min_pvalue_list_filtered, "results/acLDL/eQTLs/acLDL_rasqual_min_pvalues_filtered.rds")


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


### Experiment with multiple indepedent QTLs per gene
#Identify multiple independet SNPs per gene hits
vcf_file = readRDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.rds")
nominal_hits = dplyr::filter(ctrl_eigen_pvalue, p_eigen < 0.01)
eigenMT_test_number = dplyr::select(nominal_hits, gene_id, n_tests)

gene_df = data_frame(gene_id = "ENSG00000091490")
gene_ranges = constructGeneRanges(dplyr::select(nominal_hits, gene_id), acldl_list$gene_metadata, cis_window = 5e5)
tabix_data = tabixFetchGenes(gene_ranges[1:10], "results/acLDL/rasqual/Ctrl_500kb.sorted.txt.gz")
df = plyr::ldply(tabix_data, .id = NULL) %>% tbl_df()
df = dplyr::left_join(df, eigenMT_test_number, by = "gene_id")
df = dplyr::group_by(df, gene_id, snp_id) %>% dplyr::mutate(p_eigen = p.adjust(p_nominal, method = "bonferroni", n_tests)) %>% dplyr::ungroup()
p_values = dplyr::filter(df,p_eigen < 0.01) %>% dplyr::arrange(-chisq)

d = dplyr::filter(p_values, gene_id == "ENSG00000117226")
c = dplyr::filter(p_values, gene_id == "ENSG00000178209")


