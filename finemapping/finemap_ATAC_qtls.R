library("devtools")
library("purrr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")


#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#List of ATAC tabix files
atac_tabix_list = list(naive = "../macrophage-chromatin/results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz",
                       IFNg = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz",
                       SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz",
                       IFNg_SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz")

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Fetch cis SNPs for all caQTLs
naive_pvalues = rasqualTools::constructGeneRanges(min_pvalues_hits$naive, atac_list$gene_metadata, cis_window = 5e4) %>%
  rasqualTools::tabixFetchGenes(qtl_granges, atac_tabix_list$naive)
ifng_pvalues = rasqualTools::constructGeneRanges(min_pvalues_hits$IFNg, atac_list$gene_metadata, cis_window = 5e4) %>%
  rasqualTools::tabixFetchGenes(atac_tabix_list$IFNg)
sl1344_pvalues = rasqualTools::constructGeneRanges(min_pvalues_hits$SL1344, atac_list$gene_metadata, cis_window = 5e4) %>%
  rasqualTools::tabixFetchGenes(atac_tabix_list$SL1344)
ifng_sl1344_pvalues = rasqualTools::constructGeneRanges(min_pvalues_hits$IFNg_SL1344, atac_list$gene_metadata, cis_window = 5e4) %>%
  rasqualTools::tabixFetchGenes(atac_tabix_list$IFNg_SL1344)

#Compile into a list
pvalue_list = list(naive = naive_pvalues, IFNg = ifng_pvalues, SL1344 = sl1344_pvalues, IFNg_SL1344 = ifng_sl1344_pvalues)
saveRDS(pvalue_list, "results/ATAC/QTLs/rasqual_QTL_pvalues.rds")
pvalue_list = readRDS("results/ATAC/QTLs/rasqual_QTL_pvalues.rds")

printGene <- function(x) {
  print(x$gene_id[1])
  return(x)
  }

#Construct credible sets
naive_cs = purrr::map(pvalue_list$naive, ~dplyr::arrange(.,p_nominal) %>%
                       printGene() %>%
                       addR2FromLead(vcf_file$genotypes) %>% 
                       dplyr::filter(R2 > 0.8))
saveRDS(naive_cs, "results/ATAC/QTLs/rasqual_naive_cs.rds")


ifng_cs = purrr::map(pvalue_list$IFNg, ~dplyr::arrange(.,p_nominal) %>%
                         printGene() %>%
                         addR2FromLead(vcf_file$genotypes) %>% 
                         dplyr::filter(R2 > 0.8))
saveRDS(ifng_cs, "results/ATAC/QTLs/rasqual_IFNg_cs.rds")


sl1344_cs = purrr::map(pvalue_list$SL1344, ~dplyr::arrange(.,p_nominal) %>%
                              printGene() %>%
                              addR2FromLead(vcf_file$genotypes) %>% 
                              dplyr::filter(R2 > 0.8))
saveRDS(sl1344_cs, "results/ATAC/QTLs/rasqual_SL1344_cs.rds")


ifng_sl1344_cs = purrr::map(pvalue_list$IFNg_SL1344, ~dplyr::arrange(.,p_nominal) %>%
                       printGene() %>%
  addR2FromLead(vcf_file$genotypes) %>% 
  dplyr::filter(R2 > 0.8))
saveRDS(ifng_sl1344_cs, "results/ATAC/QTLs/rasqual_IFNg_SL1344_cs.rds")


#Save credible sets to disk
credible_set_list = list(naive = naive_cs, IFNg = ifng_cs, SL1344 = sl1344_cs, IFNg_SL1344 = ifng_sl1344_cs)
saveRDS(credible_set_list, "results/ATAC/QTLs/rasqual_credible_sets.rds")
cs = readRDS("results/ATAC/QTLs/rasqual_credible_sets.rds")
