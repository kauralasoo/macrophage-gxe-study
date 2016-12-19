library("plyr")
library("dplyr")
library("devtools")
library("Rsamtools")
library("purrr")
load_all("../seqUtils/")

##### Import RASQUAL results #####
#Nominal RASQUAL p-values
naive_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/SL1344/rasqual/naive_500kb.eigenMT.txt")
IFNg_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/SL1344/rasqual/IFNg_500kb.eigenMT.txt")
SL1344_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/SL1344/rasqual/SL1344_500kb.eigenMT.txt")
IFNg_SL1344_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/SL1344/rasqual/IFNg_SL1344_500kb.eigenMT.txt")

min_pvalue_list = list(naive = naive_eigenMT, IFNg = IFNg_eigenMT, 
                       SL1344 = SL1344_eigenMT, IFNg_SL1344 = IFNg_SL1344_eigenMT)

#Permutation Rasqual p-values
naive_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/SL1344/rasqual/random_permutation/naive_500kb.eigenMT.txt")
IFNg_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/SL1344/rasqual/random_permutation/IFNg_500kb.eigenMT.txt")
SL1344_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/SL1344/rasqual/random_permutation/SL1344_500kb.eigenMT.txt")
IFNg_SL1344_eigenMT = eigenMTImportResults("/Volumes/JetDrive/databases/SL1344/rasqual/random_permutation/IFNg_SL1344_500kb.eigenMT.txt")

min_pvalue_list_random = list(naive = naive_eigenMT, IFNg = IFNg_eigenMT, 
                       SL1344 = SL1344_eigenMT, IFNg_SL1344 = IFNg_SL1344_eigenMT)

#Calculate empirical FDR threshold
rasqual_fdr_thres_list = purrr::map2(min_pvalue_list, min_pvalue_list_random, 
                                     ~dplyr::mutate(.x, fdr_thresh = getFDR(p_eigen, .y$p_eigen, alpha = 0.1)))
saveRDS(rasqual_fdr_thres_list, "results/SL1344/eQTLs/rasqual_min_pvalues.rds")

###### FastQTL ######
#Extract the number of tests performed by eigenMT for each peak
n_tests = map(min_pvalue_list, ~dplyr::select(.,gene_id, n_tests)) %>% reduce(rbind) %>% unique()

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import nominal p-values
naive_qtls = importFastQTLTable("/Volumes/JetDrive/databases/SL1344/fastqtl/naive_500kb_permuted.txt.gz") %>% 
  enrichFastQTLPvalues(gene_name_map) %>% fastQTLCorrectEigenMT(n_tests)
IFNg_qtls = importFastQTLTable("/Volumes/JetDrive/databases/SL1344/fastqtl/IFNg_500kb_permuted.txt.gz") %>% 
  enrichFastQTLPvalues(gene_name_map) %>% fastQTLCorrectEigenMT(n_tests)
SL1344_qtls = importFastQTLTable("/Volumes/JetDrive/databases/SL1344/fastqtl/SL1344_500kb_permuted.txt.gz") %>% 
  enrichFastQTLPvalues(gene_name_map) %>% fastQTLCorrectEigenMT(n_tests)
IFNg_SL1344_qtls = importFastQTLTable("/Volumes/JetDrive/databases/SL1344/fastqtl/IFNg_SL1344_500kb_permuted.txt.gz") %>% 
  enrichFastQTLPvalues(gene_name_map) %>% fastQTLCorrectEigenMT(n_tests)

qtl_list = list(naive = naive_qtls, IFNg = IFNg_qtls, SL1344 = SL1344_qtls, IFNg_SL1344 = IFNg_SL1344_qtls)

#Import permutation p-values
naive_qtls = importFastQTLTable("/Volumes/JetDrive/databases/SL1344/fastqtl/random_permutation/naive_500kb_permuted.txt.gz") %>% 
  enrichFastQTLPvalues(gene_name_map) %>% fastQTLCorrectEigenMT(n_tests)
IFNg_qtls = importFastQTLTable("/Volumes/JetDrive/databases/SL1344/fastqtl/random_permutation/IFNg_500kb_permuted.txt.gz") %>% 
  enrichFastQTLPvalues(gene_name_map) %>% fastQTLCorrectEigenMT(n_tests)
SL1344_qtls = importFastQTLTable("/Volumes/JetDrive/databases/SL1344/fastqtl/random_permutation/SL1344_500kb_permuted.txt.gz") %>% 
  enrichFastQTLPvalues(gene_name_map) %>% fastQTLCorrectEigenMT(n_tests)
IFNg_SL1344_qtls = importFastQTLTable("/Volumes/JetDrive/databases/SL1344/fastqtl/random_permutation/IFNg_SL1344_500kb_permuted.txt.gz") %>% 
  enrichFastQTLPvalues(gene_name_map) %>% fastQTLCorrectEigenMT(n_tests)

qtl_list_random = list(naive = naive_qtls, IFNg = IFNg_qtls, SL1344 = SL1344_qtls, IFNg_SL1344 = IFNg_SL1344_qtls)

#Calculate empirical FDR threshold in each condition
fastqtl_fdr_thres_list = purrr::map2(qtl_list, qtl_list_random, 
                                     ~dplyr::mutate(.x, fdr_thresh = getFDR(p_eigen, .y$p_eigen)))
saveRDS(fastqtl_fdr_thres_list, "results/SL1344/eQTLs/fastqtl_min_pvalues.rds")


