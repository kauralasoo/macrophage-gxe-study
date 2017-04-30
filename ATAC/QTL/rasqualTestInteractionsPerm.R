library("devtools")
library("dplyr")
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
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
gene_name_map = dplyr::select(atac_list$gene_metadata, gene_id, gene_name)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Load p-values from disk
rasqual_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

#Permute conditions within individuals
perm_conditions = dplyr::group_by(atac_list$sample_metadata, donor) %>% 
  dplyr::mutate(condition_new = sample(condition_char)) %>% 
  dplyr::select(donor, condition_char, condition_new) %>% dplyr::ungroup() %>% 
  dplyr::mutate(perm_sample_id = paste(donor, condition_new, "ATAC", sep = "_"))

#Replace column names of the cqn matrix
cqn_perm = atac_list$cqn
colnames(cqn_perm) = perm_conditions$perm_sample_id

#Naive vs IFNg
covariate_names = c("sex_binary", "cqn_PC1", "cqn_PC2", "cqn_PC3")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions (lm)
interaction_results_lm = testMultipleInteractions(tbl_df(filtered_pairs), cqn_perm, 
                                                  atac_list$sample_metadata, 
                                                  vcf_file, formula_qtl, formula_interaction, id_field_separator = "-")
interaction_df_lm = postProcessInteractionPvalues(interaction_results_lm, id_field_separator = "-")

#Save to disk
saveRDS(interaction_df_lm, paste0("results/ATAC/QTLs/perm/interaction_results_lm_", perm_number,".rds"))

