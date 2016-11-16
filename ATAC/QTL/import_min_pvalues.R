library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("Rsamtools")
library("purrr")

#Import data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Import SNP coordinates
snp_coords = importVariantInformation("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

### Genotypes optimised by RASQUAL
#Extract minimal p-value for each condition
naive_eigen_pvalue = eigenMTImportResults("/Volumes/JetDrive/databases/ATAC/rasqual/naive_50kb.eigenMT.txt.gz")
ifng_eigen_pvalue = eigenMTImportResults("/Volumes/JetDrive/databases/ATAC/rasqual/IFNg_50kb.eigenMT.txt.gz")
sl1344_eigen_pvalue = eigenMTImportResults("/Volumes/JetDrive/databases/ATAC/rasqual/SL1344_50kb.eigenMT.txt.gz")
ifng_sl1344_eigen_pvalue = eigenMTImportResults("/Volumes/JetDrive/databases/ATAC/rasqual/IFNg_SL1344_50kb.eigenMT.txt.gz")

min_pvalue_list = list(naive = naive_eigen_pvalue,
                       IFNg = ifng_eigen_pvalue,
                       SL1344 = sl1344_eigen_pvalue,
                       IFNg_SL1344 = ifng_sl1344_eigen_pvalue)
saveRDS(min_pvalue_list, "results/ATAC/QTLs/rasqual_min_pvalues.rds")

#Extract the number of tests performed by eigenMT for each peak
n_tests = map(min_pvalue_list, ~dplyr::select(.,gene_id, n_tests)) %>% reduce(rbind) %>% unique()

#Find minimal p-values from fastQTL results
naive_fqtl = importFastQTLTable("results/ATAC/fastqtl/output/naive_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)
ifng_fqtl = importFastQTLTable("results/ATAC/fastqtl/output/IFNg_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)
sl1344_fqtl = importFastQTLTable("results/ATAC/fastqtl/output/SL1344_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)
ifng_sl1344_fqtl = importFastQTLTable("results/ATAC/fastqtl/output/IFNg_SL1344_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)

fastqtl_pvalue_list = list(naive = naive_fqtl,
                           IFNg = ifng_fqtl,
                           SL1344 = sl1344_fqtl, 
                           IFNg_SL1344 = ifng_sl1344_fqtl)
saveRDS(fastqtl_pvalue_list, "results/ATAC/QTLs/fastqtl_min_pvalues.rds")


#Import min p-values from fastqtl 100kb and 50kb windows centred around the peak
fastqtl_100kb_list = list(
  naive = "/Volumes/JetDrive/databases/ATAC/fastqtl/peak_centre/naive_100kb_cqn_perm.txt.gz",
  IFNg = "/Volumes/JetDrive/databases/ATAC/fastqtl/peak_centre/IFNg_100kb_cqn_perm.txt.gz",
  SL1344 = "/Volumes/JetDrive/databases/ATAC/fastqtl/peak_centre/SL1344_100kb_cqn_perm.txt.gz",
  IFNg_SL1344 = "/Volumes/JetDrive/databases/ATAC/fastqtl/peak_centre/IFNg_SL1344_100kb_cqn_perm.txt.gz")
fastqtl_100kb_min_pvalues = purrr::map(fastqtl_100kb_list, ~importFastQTLTable(.))
saveRDS(fastqtl_100kb_min_pvalues, "results/ATAC/QTLs/fastqtl_min_pvalues_100kb.rds")

fastqtl_50kb_list = list(
  naive = "/Volumes/JetDrive/databases/ATAC/fastqtl/peak_centre/naive_50kb_cqn_perm.txt.gz",
  IFNg = "/Volumes/JetDrive/databases/ATAC/fastqtl/peak_centre/IFNg_50kb_cqn_perm.txt.gz",
  SL1344 = "/Volumes/JetDrive/databases/ATAC/fastqtl/peak_centre/SL1344_50kb_cqn_perm.txt.gz",
  IFNg_SL1344 = "/Volumes/JetDrive/databases/ATAC/fastqtl/peak_centre/IFNg_SL1344_50kb_cqn_perm.txt.gz")
fastqtl_50kb_min_pvalues = purrr::map(fastqtl_50kb_list, ~importFastQTLTable(.))
saveRDS(fastqtl_50kb_min_pvalues, "results/ATAC/QTLs/fastqtl_min_pvalues_50kb.rds")



#TODO: Move this to replicability analysis
#Calculate Pi1
pi1_stat = calculatePairwisePi1(fastqtl_pvalue_list)
write.table(pi1_stat, "results/ATAC/QTLs/properties/fastqtl_pi1_results.txt", sep = "\t", quote = FALSE)
pi1_stat_tidy = calculatePairwisePi1(fastqtl_pvalue_list, tidy = TRUE)
write.table(pi1_stat_tidy, "results/ATAC/QTLs/properties/fastqtl_pi1_results_tidy.txt", sep = "\t", quote = FALSE)



