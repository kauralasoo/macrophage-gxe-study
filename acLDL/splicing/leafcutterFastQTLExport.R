library("readr")
library("dplyr")
library("tidyr")
library("limma")
library("purrr")
load_all("../seqUtils/")

#Import proportion data
prop_list = readRDS("results/acLDL/acLDL_combined_proportions.row_quantile.rds")

#Put each condition into a separate list
condition_list = list(Ctrl = "Ctrl",AcLDL = "AcLDL")
prop_conditions = purrr::map(condition_list, ~extractConditionFromExpressionList(.,prop_list))
prop_conditions_renamed = lapply(prop_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#### Export data for FastQTL ####
fastqtl_genepos = constructFastQTLGenePos(prop_conditions_renamed$Ctrl$gene_metadata)
norm_prop_list = lapply(prop_conditions_renamed, function(x){x$cqn})
fastqtl_norm_prop_list = lapply(norm_prop_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_norm_prop_list, "results/acLDL/leafcutter/fastqtl_input/", file_suffix = "norm_prop")

#Save covariates
covariate_names = c("genotype_id", "norm_PC1", "norm_PC2", "norm_PC3", "norm_PC4","norm_PC5","sex_binary")
covariate_list = lapply(prop_conditions_renamed, function(x, names){x$sample_metadata[,names]}, covariate_names)
fastqtl_covariates = lapply(covariate_list, fastqtlMetadataToCovariates)
saveFastqtlMatrices(fastqtl_covariates, "results/acLDL/leafcutter/fastqtl_input/", file_suffix = "covariates_prop")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:250), n = 250)
write.table(chunks_matrix, "results/acLDL/leafcutter/fastqtl_input/all_chunk_table.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")
