library("purrr")
library("dplyr")
library("ggplot2")

importMinPvaluePerGene <- function(perm_number, type = "lm"){
  if(type == "lm"){
    path = paste0("results/SL1344/eQTLs/perm/interaction_results_lm_", perm_number, ".rds")
  } else {
    path = paste0("results/SL1344/eQTLs/perm/interaction_results_lme4_", perm_number, ".rds")
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
interaction_df_lm = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues.rds") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene_id, snp_id, p_nominal)
interaction_df_lme4 = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues_lme4.rds") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene_id, snp_id, p_nominal)

#Construct a list of permutation numbers
perm_numbers = as.list(c(1:1000))
names(perm_numbers) = c(1:1000)

#Import permutation p-values
perm_pvalues_lm = purrr::map_df(perm_numbers, ~importMinPvaluePerGene(.,type = "lm"), .id = "perm_number")
saveRDS(perm_pvalues_lm, "results/SL1344/eQTLs/SL1344_interaction_pvalues.permuted.rds")
perm_pvalues_lm = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues.permuted.rds")

#Colculate empirical p-values for each gene
empirical_pvalues_lm = dplyr::rename(perm_pvalues_lm, p_all = p_nominal) %>%
  dplyr::left_join(interaction_df_lm, by = "gene_id") %>%
  dplyr::mutate(is_larger = as.numeric(p_nominal > p_all)) %>%
  dplyr::group_by(gene_id, snp_id) %>%
  dplyr::summarise(p_perm = (sum(is_larger) + 1)/1000) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_perm = pmin(p_perm, 1)) %>%
  dplyr::mutate(p_fdr = p.adjust(p_perm, method = "fdr"))
saveRDS(empirical_pvalues_lm, "results/SL1344/eQTLs/SL1344_interaction_pvalues_lm.empirical.rds")


#Import permutation p-values from lme4
perm_pvalues_lme4 = purrr::map_df(perm_numbers, ~importMinPvaluePerGene(.,type = "lme4"), .id = "perm_number")
saveRDS(perm_pvalues_lme4, "results/SL1344/eQTLs/SL1344_interaction_pvalues_lme4.permuted.rds")

empirical_pvalues_lme4 = dplyr::rename(perm_pvalues_lme4, p_all = p_nominal) %>%
  dplyr::left_join(interaction_df_lme4, by = "gene_id") %>%
  dplyr::mutate(is_larger = as.numeric(p_nominal > p_all)) %>%
  dplyr::group_by(gene_id, snp_id) %>%
  dplyr::summarise(p_perm = (sum(is_larger) + 1)/1000) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_perm = pmin(p_perm, 1)) %>%
  dplyr::mutate(p_fdr = p.adjust(p_perm, method = "fdr"))
saveRDS(empirical_pvalues_lme4, "results/SL1344/eQTLs/SL1344_interaction_pvalues_lme4.empirical.rds")



