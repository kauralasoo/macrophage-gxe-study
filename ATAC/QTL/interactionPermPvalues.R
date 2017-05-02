library("purrr")
library("dplyr")
library("ggplot2")
library("tidyr")

importMinPvaluePerGene <- function(perm_number, type = "lm"){
  print(perm_number)
  if(type == "lm"){
    path = paste0("results/ATAC/QTLs/perm/interaction_results_lm_", perm_number, ".rds")
  } else {
    path = paste0("results/ATAC/QTLs/perm/interaction_results_lme4_", perm_number, ".rds")
  }
  rds = readRDS(path) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::arrange(gene_id, p_nominal) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene_id, p_nominal)
  return(rds)
}

#Import p-values from the nominal test
interaction_df_lm = readRDS("results/ATAC/QTLs/rasqual_interaction_results.rds") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene_id, snp_id, p_nominal)

#Construct a list of permutation numbers
perm_files = data_frame(file_name = list.files("results/ATAC/QTLs/perm"))
perm_number_df = tidyr::separate(perm_files, file_name, c("prefix", "suffix"), sep = "_lm_") %>%
  tidyr::separate(suffix, c("perm_number", "suffix"), sep = ".rds")
perm_numbers = perm_number_df$perm_number
names(perm_numbers) = perm_numbers

#Import permutation p-values
perm_pvalues_lm = purrr::map_df(perm_numbers, ~importMinPvaluePerGene(.,type = "lm"), .id = "perm_number")
saveRDS(perm_pvalues_lm, "results/ATAC/QTLs/ATAC_interaction_pvalues.permuted.rds")

#Colculate empirical p-values for each gene
empirical_pvalues_lm = dplyr::rename(perm_pvalues_lm, p_all = p_nominal) %>%
  dplyr::left_join(interaction_df_lm, by = "gene_id") %>%
  dplyr::mutate(is_larger = as.numeric(p_nominal > p_all)) %>%
  dplyr::group_by(gene_id, snp_id) %>%
  dplyr::summarise(p_perm = (sum(is_larger) + 1)/1000) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_perm = pmin(p_perm, 1)) %>%
  dplyr::mutate(p_fdr = p.adjust(p_perm, method = "fdr"))
saveRDS(empirical_pvalues_lm, "results/ATAC/QTLs/ATAC_interaction_pvalues_lm.empirical.rds")


