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


#Import fastQTL p-values from disk
ctrl_perm_pvalue = importFastQTLTable("results/acLDL/fastqtl/output/Ctrl_permuted.txt.gz")
acldl_perm_pvalue = importFastQTLTable("results/acLDL/fastqtl/output/AcLDL_permuted.txt.gz")
perm_pvalue_list = list(Ctrl = ctrl_perm_pvalue, AcLDL = acldl_perm_pvalue)
saveRDS(perm_pvalue_list, "results/acLDL/eQTLs/fastqtl_min_pvalues.rds")
a = readRDS("results/acLDL/eQTLs/fastqtl_min_pvalues.rds")

#Import QTLtools p-values from disk
ctrl_perm = importQTLtoolsTable("processed/acLDL/fastqtl_output/featureCounts/Ctrl.permuted.txt.gz")
acldl_perm = importQTLtoolsTable("processed/acLDL/fastqtl_output/featureCounts/AcLDL.permuted.txt.gz")
perm_pvalue_list = list(Ctrl = ctrl_perm, AcLDL = acldl_perm)
saveRDS(perm_pvalue_list, "results/acLDL/eQTLs/qtltools_min_pvalues.rds")

#Compare the correlation betweeb rasqual and fastQTL p-values
fastqtl_vs_rasqual = ggplot(ctrl_matched, aes(x = -log(p_linear, 10), y = -log(p_rasqual,10))) + geom_point()
ggsave("results/acLDL/eQTLs/fastqtl_vs_rasqual.pdf", plot = fastqtl_vs_rasqual)

#Export lead-qQTLs
dplyr::left_join(ctrl_eigen_pvalue, gene_name_map, by = "gene_id") %>%
  write.table("results/acLDL/eQTLs/ctrl_rasqual_hits.txt", quote = FALSE, row.names = FALSE, sep = "\t")
dplyr::left_join(acldl_eigen_pvalue, gene_name_map, by = "gene_id") %>%
  write.table("results/acLDL/eQTLs/acLDL_rasqual_hits.txt", quote = FALSE, row.name = FALSE, sep = "\t")

