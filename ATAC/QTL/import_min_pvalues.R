library("plyr")
library("dplyr")
library("devtools")
library("Rsamtools")
library("purrr")
load_all("../seqUtils/")

#Import data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Import SNP coordinates
snp_coords = importVariantInformation("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

### RASQUAL nominal run
#Extract minimal p-value for each condition
naive_eigen_pvalue = eigenMTImportResults("~/databases/ATAC/rasqual/naive_50kb.eigenMT.txt.gz")
ifng_eigen_pvalue = eigenMTImportResults("~/databases/ATAC/rasqual/IFNg_50kb.eigenMT.txt.gz")
sl1344_eigen_pvalue = eigenMTImportResults("~/databases/ATAC/rasqual/SL1344_50kb.eigenMT.txt.gz")
ifng_sl1344_eigen_pvalue = eigenMTImportResults("~/databases/ATAC/rasqual/IFNg_SL1344_50kb.eigenMT.txt.gz")

min_pvalue_list = list(naive = naive_eigen_pvalue,
                       IFNg = ifng_eigen_pvalue,
                       SL1344 = sl1344_eigen_pvalue,
                       IFNg_SL1344 = ifng_sl1344_eigen_pvalue)

### RASQUAL permutation run
#Extract minimal p-value for each condition
naive_eigen_pvalue = eigenMTImportResults("~/databases/ATAC/rasqual/random_permutation/naive_50kb.eigenMT.txt.gz")
ifng_eigen_pvalue = eigenMTImportResults("~/databases/ATAC/rasqual/random_permutation/IFNg_50kb.eigenMT.txt.gz")
sl1344_eigen_pvalue = eigenMTImportResults("~/databases/ATAC/rasqual/random_permutation/SL1344_50kb.eigenMT.txt.gz")
ifng_sl1344_eigen_pvalue = eigenMTImportResults("~/databases/ATAC/rasqual/random_permutation/IFNg_SL1344_50kb.eigenMT.txt.gz")

min_pvalue_list_random = list(naive = naive_eigen_pvalue,
                       IFNg = ifng_eigen_pvalue,
                       SL1344 = sl1344_eigen_pvalue,
                       IFNg_SL1344 = ifng_sl1344_eigen_pvalue)

#Calculate empirical FDR threshold
rasqual_fdr_thres_list = purrr::map2(min_pvalue_list, min_pvalue_list_random, 
                                     ~dplyr::mutate(.x, fdr_thresh = getFDR(p_eigen, .y$p_eigen)))
saveRDS(rasqual_fdr_thres_list, "results/ATAC/QTLs/rasqual_min_pvalues.rds")
purrr::map(rasqual_fdr_thres_list, ~dplyr::filter(.,p_eigen < fdr_thresh))

#### FastQTL ####

#Extract the number of tests performed by eigenMT for each peak
n_tests = map(min_pvalue_list, ~dplyr::select(.,gene_id, n_tests)) %>% reduce(rbind) %>% unique()

#Find minimal p-values from fastQTL results
naive_fqtl = importFastQTLTable("~/databases/ATAC/fastqtl/naive_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)
ifng_fqtl = importFastQTLTable("~/databases/ATAC/fastqtl/IFNg_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)
sl1344_fqtl = importFastQTLTable("~/databases/ATAC/fastqtl/SL1344_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)
ifng_sl1344_fqtl = importFastQTLTable("~/databases/ATAC/fastqtl/IFNg_SL1344_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)

#Merge to a list
fastqtl_pvalue_list = list(naive = naive_fqtl,
                           IFNg = ifng_fqtl,
                           SL1344 = sl1344_fqtl, 
                           IFNg_SL1344 = ifng_sl1344_fqtl)

#Import min p-values from the permutation run
naive_fqtl = importFastQTLTable("~/databases/ATAC/fastqtl/random_permutation/naive_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)
ifng_fqtl = importFastQTLTable("~/databases/ATAC/fastqtl/random_permutation/IFNg_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)
sl1344_fqtl = importFastQTLTable("~/databases/ATAC/fastqtl/random_permutation/SL1344_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)
ifng_sl1344_fqtl = importFastQTLTable("~/databases/ATAC/fastqtl/random_permutation/IFNg_SL1344_50kb_cqn_perm.txt.gz") %>% 
  fastQTLCorrectEigenMT(n_tests)

#Merge to a list
fastqtl_pvalue_list_random = list(naive = naive_fqtl,
                           IFNg = ifng_fqtl,
                           SL1344 = sl1344_fqtl, 
                           IFNg_SL1344 = ifng_sl1344_fqtl)
saveRDS(fastqtl_pvalue_list_random, "results/ATAC/QTLs/fastqtl_random_pvalues.rds")


#Calculate empirical FDR threshold in each condition
fastqtl_fdr_thres_list = purrr::map2(fastqtl_pvalue_list, fastqtl_pvalue_list_random, 
                        ~dplyr::mutate(.x, fdr_thresh = getFDR(p_eigen, .y$p_eigen)))
saveRDS(fastqtl_fdr_thres_list, "results/ATAC/QTLs/fastqtl_min_pvalues.rds")

purrr::map(fastqtl_fdr_thres_list, ~dplyr::filter(.,p_eigen < fdr_thresh))

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


##### Export min QTL lists for Zenodo #####
#RASQUAL
rasqual_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
write.table(rasqual_pvalues$naive, "figures/tables/caQTLs/caQTL_naive_RASQUAL_lead_only.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rasqual_pvalues$IFNg, "figures/tables/caQTLs/caQTL_IFNg_RASQUAL_lead_only.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rasqual_pvalues$SL1344, "figures/tables/caQTLs/caQTL_SL1344_RASQUAL_lead_only.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rasqual_pvalues$IFNg_SL1344, "figures/tables/caQTLs/caQTL_IFNg_SL1344_RASQUAL_lead_only.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#FastQTL
fastqtl_pvalues = readRDS("results/ATAC/QTLs/fastqtl_min_pvalues.rds")
write.table(fastqtl_pvalues$naive, "figures/tables/caQTLs/caQTL_naive_FastQTL_lead_only.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fastqtl_pvalues$IFNg, "figures/tables/caQTLs/caQTL_IFNg_FastQTL_lead_only.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fastqtl_pvalues$SL1344, "figures/tables/caQTLs/caQTL_SL1344_FastQTL_lead_only.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fastqtl_pvalues$IFNg_SL1344, "figures/tables/caQTLs/caQTL_IFNg_SL1344_FastQTL_lead_only.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

