library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("Rsamtools")
library("purrr")

#Import data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")

#Import SNP coordinates
snp_coords = read.table("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.snp_coords.txt", stringsAsFactors = FALSE)
colnames(snp_coords) = c("chr", "pos", "snp_id")

### Genotypes optimised by RASQUAL
#Extract minimal p-value for each condition
naive_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output/naive_100kb/naive_50kb.eigenMT.txt.gz")
ifng_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output/IFNg_100kb/IFNg_50kb.eigenMT.txt.gz")
sl1344_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output/SL1344_100kb/SL1344_50kb.eigenMT.txt.gz")
ifng_sl1344_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_50kb.eigenMT.txt.gz")

min_pvalue_list = list(naive = naive_eigen_pvalue,
                       IFNg = ifng_eigen_pvalue,
                       SL1344 = sl1344_eigen_pvalue,
                       IFNg_SL1344 = ifng_sl1344_eigen_pvalue)
saveRDS(min_pvalue_list, "results/ATAC/QTLs/rasqual_min_pvalues.rds")

#Fixed genotypes
#Extract minimal p-value for each condition
naive_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output_fixed_gt//naive_100kb/naive_50kb.eigenMT.txt.gz")
ifng_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output_fixed_gt/IFNg_100kb/IFNg_50kb.eigenMT.txt.gz")
sl1344_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output_fixed_gt/SL1344_100kb/SL1344_50kb.eigenMT.txt.gz")
ifng_sl1344_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output_fixed_gt/IFNg_SL1344_100kb/IFNg_SL1344_50kb.eigenMT.txt.gz")

min_pvalue_list = list(naive = naive_eigen_pvalue,
                       IFNg = ifng_eigen_pvalue,
                       SL1344 = sl1344_eigen_pvalue,
                       IFNg_SL1344 = ifng_sl1344_eigen_pvalue)
saveRDS(min_pvalue_list, "results/ATAC/QTLs/rasqual_min_pvalues.fixed_gt.rds")


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

#Calculate Pi1
pi1_stat = calculatePairwisePi1(fastqtl_pvalue_list)
write.table(pi1_stat, "results/ATAC/QTLs/fastqtl_pi1_results.txt", sep = "\t", quote = FALSE)
