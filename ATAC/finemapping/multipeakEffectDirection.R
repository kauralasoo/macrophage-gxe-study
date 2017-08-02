library("devtools")
library("dplyr")
library("purrr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

#Helper functions
fetchRasqualSNPs <- function(snp_ids, snpspos, summary_list){
  #Construct GRanges object of SNP positions
  selected_snps = dplyr::filter(snpspos, snpid %in% snp_ids) %>%
    dplyr::transmute(snp_id = snpid, seqnames = chr, start = pos, end = pos, strand = "*") %>%
    dataFrameToGRanges()
  
  selected_pvalues = lapply(summary_list, function(tabix, snps) rasqualTools::tabixFetchSNPs(snps, tabix), selected_snps)
  return(selected_pvalues)
}

addEffectSizes <- function(master_dependent_pairs, rasqual_SNPs){
  master_effects = dplyr::select(master_dependent_pairs, master_id, snp_id) %>% 
    unique() %>% 
    dplyr::left_join(rasqual_SNPs, by = c("master_id" = "gene_id", "snp_id")) %>% 
    dplyr::transmute(master_id, snp_id, master_beta = beta, master_p = p_nominal)
  
  dependent_effects = dplyr::select(master_dependent_pairs, dependent_id, snp_id, master_id) %>% 
    unique() %>% 
    dplyr::left_join(rasqual_SNPs, by = c("dependent_id" = "gene_id", "snp_id")) %>% 
    dplyr::transmute(dependent_id, master_id, snp_id, dependent_beta = beta, dependent_p = p_nominal)
  
  effects_df = dplyr::left_join(dependent_effects, master_effects, by = c("master_id", "snp_id")) %>%
    dplyr::mutate(sign_diff = abs(sign(master_beta) - sign(dependent_beta)))
  return(effects_df)
}



#Import atac data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import Master-dependent pairs
result_list = readRDS("results/ATAC/QTLs/qtl_peak_type_assignment.rds")
master_dependent_pairs = result_list$dependents$unique_masters

#Fetch rasqual SNPs from disk
atac_selected_pvalues = fetchRasqualSNPs(unique(master_dependent_pairs$snp_id), vcf_file$snpspos, qtlResults()$atac_rasqual)

#Add effect sizes
effects_df_list = purrr::map(atac_selected_pvalues, ~addEffectSizes(master_dependent_pairs, .) %>%
                               dplyr::filter(master_p < 1e-4, dependent_p < 1e-4))
purrr::map_df(effects_df_list, ~ table(.$sign_diff)[2]/sum(table(.$sign_diff)))

