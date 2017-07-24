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

