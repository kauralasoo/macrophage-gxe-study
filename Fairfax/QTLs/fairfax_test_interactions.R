library("dplyr")
library("tidyr")
library("purrr")
library("UpSetR")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import SummarizedExperiment
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")

#keep only sample that are present in all conditions
selected_samples = (colData(se_fairfax) %>% tbl_df2() %>% dplyr::filter(present_in_all))$sample_id
se_shared = se_fairfax[,selected_samples]

#Extract expression and metadata
expression_mat = assays(se_shared)$exprs
sample_metadata = colData(se_shared) %>% tbl_df2() %>% 
  dplyr::rename(donor = donor_id)
gene_metadata = rowData(se_shared) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id)

#Calculate covariates (PCs)
#Extract conditions
condition_list = idVectorToList(c("CD14","IFN","LPS2","LPS24"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,se_shared))

#Calculate PCs in each condition
sample_meta = purrr::map(event_conditions, ~tbl_df2(colData(.)))
exprs_list = purrr::map(event_conditions, ~assays(.)$exprs)
covariates_list = purrr::map2(exprs_list, sample_meta, 
                              ~performPCA(.x, .y, n_pcs = 6, feature_id = "sample_id")$pca_matrix %>%
                                dplyr::select(sample_id, PC1, PC2, PC3, PC4, PC5, PC6))
sample_meta_cov = dplyr::left_join(sample_metadata, purrr::map_df(covariates_list, identity), by = "sample_id")


#Import QTLs
qtl_min_pvalues = readRDS("results/Fairfax/fairfax_qtl_min_pvalues.rds")

#Identify all qtl pairs
min_pvalue_hits = purrr::map(qtl_min_pvalues$shared, ~dplyr::filter(.,p_fdr < 0.1))
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(phenotype_id, p_nominal)
joint_pairs = dplyr::transmute(min_pvalues_df, gene_id = phenotype_id, snp_id, pheno_chr) %>% unique()

#Export QTL pairs for aFC calculation
afc_pairs = dplyr::transmute(min_pvalues_df, pid = phenotype_id, sid = snp_id, sid_chr = snp_chr, sid_pos = snp_start) %>% unique()
write.table(afc_pairs, "processed/Fairfax/qtltools/input/shared/qtl_pairs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#Use a paired design to test for interaction
covariate_names = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions chr by chr
chr_list = c(1:22) %>% as.character() %>% idVectorToList() 
interaction_res_list = purrr::map(chr_list, ~testInteractionByChromosome(chr = ., joint_pairs, trait_matrix = expression_mat, 
                             sample_metadata = sample_meta_cov, 
                             qtl_formula = formula_qtl, 
                             interaction_formula = formula_interaction, 
                             id_field_separator = "-", lme4 = TRUE))
interaction_df = purrr::flatten(interaction_res_list) %>% postProcessInteractionPvalues(id_field_separator = "-")
saveRDS(interaction_df, "results/Fairfax/shared_interactions_lme4.rds")



#Perform interaction test on a random subset of 84 individuals
set.seed(42)
shared_samples = (colData(se_fairfax) %>% tbl_df2() %>% dplyr::filter(present_in_all))
shared_donors = unique(shared_samples$donor_id)
selected_samples = dplyr::filter(shared_samples, donor_id %in% sample(shared_donors, 84))
se_shared = se_fairfax[,selected_samples$sample_id]

#Extract expression and metadata
expression_mat = assays(se_shared)$exprs
sample_metadata = colData(se_shared) %>% tbl_df2() %>% 
  dplyr::rename(donor = donor_id)
gene_metadata = rowData(se_shared) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id)

#Calculate covariates (PCs)
#Extract conditions
condition_list = idVectorToList(c("CD14","IFN","LPS2","LPS24"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,se_shared))

#Calculate PCs in each condition
sample_meta = purrr::map(event_conditions, ~tbl_df2(colData(.)))
exprs_list = purrr::map(event_conditions, ~assays(.)$exprs)
covariates_list = purrr::map2(exprs_list, sample_meta, 
                              ~performPCA(.x, .y, n_pcs = 6, feature_id = "sample_id")$pca_matrix %>%
                                dplyr::select(sample_id, PC1, PC2, PC3, PC4, PC5, PC6))
sample_meta_cov = dplyr::left_join(sample_metadata, purrr::map_df(covariates_list, identity), by = "sample_id")


#Import QTLs
qtl_min_pvalues = readRDS("results/Fairfax/fairfax_qtl_min_pvalues.rds")

#Identify all qtl pairs
min_pvalue_hits = purrr::map(qtl_min_pvalues$shared_84, ~dplyr::filter(.,p_fdr < 0.1))
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(phenotype_id, p_nominal)
joint_pairs = dplyr::transmute(min_pvalues_df, gene_id = phenotype_id, snp_id, pheno_chr) %>% unique()

#Export QTL pairs for aFC calculation
afc_pairs = dplyr::transmute(min_pvalues_df, pid = phenotype_id, sid = snp_id, sid_chr = snp_chr, sid_pos = snp_start) %>% unique()
write.table(afc_pairs, "processed/Fairfax/qtltools/input/shared_84/qtl_pairs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Use a paired design to test for interaction
covariate_names = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions chr by chr
chr_list = c(1:22) %>% as.character() %>% idVectorToList() 
interaction_res_list = purrr::map(chr_list, ~testInteractionByChromosome(chr = ., joint_pairs, trait_matrix = expression_mat, 
                                                                         sample_metadata = sample_meta_cov, 
                                                                         qtl_formula = formula_qtl, 
                                                                         interaction_formula = formula_interaction, 
                                                                         id_field_separator = "-", lme4 = TRUE))
interaction_df = purrr::flatten(interaction_res_list) %>% postProcessInteractionPvalues(id_field_separator = "-")
saveRDS(interaction_df, "results/Fairfax/shared_84_interactions_lme4.rds")




#Perform interaction test on a random subset of 42 individuals
set.seed(42)
shared_samples = (colData(se_fairfax) %>% tbl_df2() %>% dplyr::filter(present_in_all))
shared_donors = unique(shared_samples$donor_id)
selected_samples = dplyr::filter(shared_samples, donor_id %in% sample(shared_donors, 42))
se_shared = se_fairfax[,selected_samples$sample_id]

#Extract expression and metadata
expression_mat = assays(se_shared)$exprs
sample_metadata = colData(se_shared) %>% tbl_df2() %>% 
  dplyr::rename(donor = donor_id)
gene_metadata = rowData(se_shared) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id)

#Calculate covariates (PCs)
#Extract conditions
condition_list = idVectorToList(c("CD14","IFN","LPS2","LPS24"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,se_shared))

#Calculate PCs in each condition
sample_meta = purrr::map(event_conditions, ~tbl_df2(colData(.)))
exprs_list = purrr::map(event_conditions, ~assays(.)$exprs)
covariates_list = purrr::map2(exprs_list, sample_meta, 
                              ~performPCA(.x, .y, n_pcs = 6, feature_id = "sample_id")$pca_matrix %>%
                                dplyr::select(sample_id, PC1, PC2, PC3, PC4, PC5, PC6))
sample_meta_cov = dplyr::left_join(sample_metadata, purrr::map_df(covariates_list, identity), by = "sample_id")


#Import QTLs
qtl_min_pvalues = readRDS("results/Fairfax/fairfax_qtl_min_pvalues.rds")

#Identify all qtl pairs
min_pvalue_hits = purrr::map(qtl_min_pvalues$shared_42, ~dplyr::filter(.,p_fdr < 0.1))
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(phenotype_id, p_nominal)
joint_pairs = dplyr::transmute(min_pvalues_df, gene_id = phenotype_id, snp_id, pheno_chr) %>% unique()

#Export QTL pairs for aFC calculation
afc_pairs = dplyr::transmute(min_pvalues_df, pid = phenotype_id, sid = snp_id, sid_chr = snp_chr, sid_pos = snp_start) %>% unique()
write.table(afc_pairs, "processed/Fairfax/qtltools/input/shared_42/qtl_pairs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Use a paired design to test for interaction
covariate_names = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions chr by chr
chr_list = c(1:22) %>% as.character() %>% idVectorToList() 
interaction_res_list = purrr::map(chr_list, ~testInteractionByChromosome(chr = ., joint_pairs, trait_matrix = expression_mat, 
                                                                         sample_metadata = sample_meta_cov, 
                                                                         qtl_formula = formula_qtl, 
                                                                         interaction_formula = formula_interaction, 
                                                                         id_field_separator = "-", lme4 = TRUE))
interaction_df = purrr::flatten(interaction_res_list) %>% postProcessInteractionPvalues(id_field_separator = "-")
saveRDS(interaction_df, "results/Fairfax/interactions/shared_42_interactions_lme4.rds")


