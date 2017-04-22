library("devtools")
library("dplyr")
library("lme4")
load_all("../seqUtils/")

batch_id = NULL

####### Get batch id from disk ######
f <- file("stdin")
open(f)
perm_number = readLines(f)
print(perm_number)
close(f)
####### END #######


#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Load p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

#Permute conditions within individuals
perm_conditions = dplyr::group_by(combined_expression_data$sample_metadata, donor) %>% 
  dplyr::mutate(condition_new = sample(condition)) %>% 
  dplyr::select(donor, condition, condition_new) %>% dplyr::ungroup() %>% 
  dplyr::mutate(perm_sample_id = paste(donor, condition_new, sep = "_"))

#Replace column names of the cqn matrix
cqn_perm = combined_expression_data$cqn
colnames(cqn_perm) = perm_conditions$perm_sample_id

#Naive vs IFNg
covariate_names = c("PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6", "sex_binary")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions (lm)
interaction_results_lm = testMultipleInteractions(tbl_df(filtered_pairs), cqn_perm, 
                                               combined_expression_data$sample_metadata, 
                                               vcf_file, formula_qtl, formula_interaction, id_field_separator = "-")
interaction_df_lm = postProcessInteractionPvalues(interaction_results_lm, id_field_separator = "-")

#Test for interactions (lme4)
##### Ctrl vs AcLDL (paired design with lme4) #####
covariate_names = c("sex_binary","PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results_lme4 = testMultipleInteractions(tbl_df(filtered_pairs), cqn_perm, 
                                               combined_expression_data$sample_metadata, 
                                               vcf_file, formula_qtl, formula_interaction, id_field_separator = "-", lme4 = TRUE)
interaction_df_lme4 = postProcessInteractionPvalues(interaction_results_lme4, id_field_separator = "-")

#Save permutatation results to disk
saveRDS(interaction_df_lm, paste0("results/acLDL/eQTLs/perm/interaction_results_lm_",perm_number,".rds"))
saveRDS(interaction_df_lme4, paste0("results/acLDL/eQTLs/perm/interaction_results_lme4_",perm_number,".rds"))


