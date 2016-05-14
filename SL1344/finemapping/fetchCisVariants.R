library("devtools")
library("purrr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")


#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#List of ATAC tabix files
rna_tabix_list = list(naive = "results/SL1344/rasqual/output/naive_500kb/naive_500kb.sorted.txt.gz",
                       IFNg = "results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.sorted.txt.gz",
                       SL1344 = "results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.sorted.txt.gz",
                       IFNg_SL1344 = "results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.sorted.txt.gz")

#Import ATAC data
rna_list = readRDS("results/SL1344/combined_expression_data_covariates.rds")
min_pvalues_list = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Fetch cis SNPs for all caQTL
pvalues_list = purrr::map2(min_pvalues_hits, rna_tabix_list, function(x,y){
  rasqualTools::constructGeneRanges(x, rna_list$gene_metadata, cis_window = 5e5) %>%
    rasqualTools::tabixFetchGenes(y)
})
saveRDS(pvalues_list, "results/SL1344/eQTLs/rasqual_QTL_pvalues.rds")

#Construct credible sets
printGene <- function(x) {
  print(x$gene_id[1])
  return(x)
}
credible_sets = purrr::map(pvalues_list, 
                           ~purrr::map(.,~dplyr::arrange(.,p_nominal) %>%
                                         printGene() %>%
                                         addR2FromLead(vcf_file$genotypes) %>% 
                                         dplyr::filter(R2 > 0.8)))
saveRDS(credible_sets, "results/SL1344/eQTLs/rasqual_credible_sets.rds")



