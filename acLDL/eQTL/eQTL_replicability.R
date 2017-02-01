library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import eQTLs
acLDL_rasqual = readRDS("results/acLDL/eQTLs/rasqual_min_pvalues.rds")
SL1344_rasqual = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")

acLDL_fastqtl = readRDS("results/acLDL/eQTLs/fastqtl_min_pvalues.rds")
SL1344_fastqtl = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")

#Take naive p-values from both studies
SL1344_shared_variants = SL1344_rasqual$naive[SL1344_rasqual$naive$snp_id %in% vcf_file$snpspos$snpid,]
acLDL_shared_genes = acLDL_rasqual$Ctrl[acLDL_rasqual$Ctrl$gene_id %in% SL1344_shared_variants$gene_id,]
SL1344_shared_variants = SL1344_shared_variants[SL1344_shared_variants$gene_id %in% acLDL_shared_genes$gene_id,]
pairwise_naive = list(Salmonella = SL1344_shared_variants, AcLDL = acLDL_shared_genes)

SL1344_shared_variants = SL1344_fastqtl$naive[SL1344_fastqtl$naive$snp_id %in% vcf_file$snpspos$snpid,]
acLDL_shared_genes = acLDL_fastqtl$Ctrl[acLDL_fastqtl$Ctrl$gene_id %in% SL1344_shared_variants$gene_id,]
SL1344_shared_variants = SL1344_shared_variants[SL1344_shared_variants$gene_id %in% acLDL_shared_genes$gene_id,]
pairwise_naive_fastqtl = list(Salmonella = SL1344_shared_variants, AcLDL = acLDL_shared_genes)

#Import the VCF file
vcf_file = readRDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.rds")

#Calculate RNA concordances
ctrl_acldl_concordance = 
  calculatePairwiseConcordance(acLDL_rasqual, ~p_nominal < 1e-8, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "RASQUAL (500kb)", class = "Ctrl vs acLDL")

naive_ctrl_concordance = 
  calculatePairwiseConcordance(pairwise_naive, ~p_nominal < 1e-8, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "RASQUAL (500kb)", class = "Ctrl vs acLDL")

ctrl_acldl_fastqtl_concordance = 
  calculatePairwiseConcordance(acLDL_rasqual, ~p_nominal < 1e-8, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (500kb)", class = "Ctrl vs acLDL")


#Calculate Pi1 stats
calculatePairwisePi1(acLDL_fastqtl, tidy = T)
calculatePairwisePi1(pairwise_naive_fastqtl, tidy = T)


                     
