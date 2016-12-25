library("devtools")
library("purrr")
library("dplyr")
library("MatrixEQTL")
library("ggplot2")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#List of ATAC tabix files
atac_tabix_list = qtlResults()$atac_rasqual

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = purrr::map(min_pvalues_list, ~dplyr::filter(., p_eigen < fdr_thresh))

#Fetch p-values for each gene
pvalues_list = purrr::map2(min_pvalues_hits, atac_tabix_list, function(x,y){
  rasqualTools::constructGeneRanges(x, atac_list$gene_metadata, cis_window = 5e4) %>%
    rasqualTools::tabixFetchGenes(y)
})
saveRDS(pvalues_list, "results/ATAC/QTLs/rasqual_QTL_pvalues.rds")
pvalues_list = readRDS("results/ATAC/QTLs/rasqual_QTL_pvalues.rds")

printGene <- function(x) {
  print(x$gene_id[1])
  return(x)
  }

#Construct credible sets
naive_cs = purrr::map(pvalues_list$naive, ~dplyr::arrange(.,p_nominal) %>%
                       printGene() %>%
                       addR2FromLead(vcf_file$genotypes) %>% 
                       dplyr::filter(R2 > 0.8))

ifng_cs = purrr::map(pvalues_list$IFNg, ~dplyr::arrange(.,p_nominal) %>%
                         printGene() %>%
                         addR2FromLead(vcf_file$genotypes) %>% 
                         dplyr::filter(R2 > 0.8))

sl1344_cs = purrr::map(pvalues_list$SL1344, ~dplyr::arrange(.,p_nominal) %>%
                              printGene() %>%
                              addR2FromLead(vcf_file$genotypes) %>% 
                              dplyr::filter(R2 > 0.8))

ifng_sl1344_cs = purrr::map(pvalues_list$IFNg_SL1344, ~dplyr::arrange(.,p_nominal) %>%
                       printGene() %>%
  addR2FromLead(vcf_file$genotypes) %>% 
  dplyr::filter(R2 > 0.8))

#Save credible sets to disk
credible_set_list = list(naive = naive_cs, IFNg = ifng_cs, SL1344 = sl1344_cs, IFNg_SL1344 = ifng_sl1344_cs)
saveRDS(credible_set_list, "results/ATAC/QTLs/rasqual_credible_sets.rds")
cs = readRDS("results/ATAC/QTLs/rasqual_credible_sets.rds")
