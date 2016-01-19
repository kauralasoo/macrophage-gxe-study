library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data.rds")

#Filter genes by expression level
mean_expression = calculateMean(combined_expression_data$cqn, as.data.frame(combined_expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_expression, 1, max) > 0))
not_Y_genes = dplyr::filter(combined_expression_data$gene_metadata, !(chromosome_name %in% c("MT","Y")))$gene_id
keep_genes = intersect(expressed_genes, not_Y_genes)
rna_expressed = extractGenesFromExpressionList(combined_expression_data, keep_genes)

#Extract separate lists for each condition
condition_names = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
rna_conditions = lapply(condition_names, extractConditionFromExpressionList, rna_expressed)

#Rename column names to genotype ids
rna_conditions_renamed = lapply(rna_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#### Set up covariates ####
#Save expression data for PEER and run PEER outside of R
cqn_list = lapply(rna_conditions_renamed, function(x){x$cqn})
savePEERData(cqn_list, "results/SL1344/PEER/input/")

#Run PEER .....

#Import PEER factors into R
cond_A_factors = importPEERFactors("results/SL1344/PEER/naive_10/factors.txt", rna_conditions_renamed$naive$sample_metadata)
cond_B_factors = importPEERFactors("results/SL1344/PEER/IFNg_10/factors.txt", rna_conditions_renamed$IFNg$sample_metadata)
cond_C_factors = importPEERFactors("results/SL1344/PEER/SL1344_10/factors.txt", rna_conditions_renamed$SL1344$sample_metadata)
cond_D_factors = importPEERFactors("results/SL1344/PEER/IFNg_SL1344_10/factors.txt", rna_conditions_renamed$IFNg_SL1344$sample_metadata)
peer_factors = rbind(cond_A_factors,cond_B_factors, cond_C_factors, cond_D_factors) %>%
  dplyr::semi_join(combined_expression_data$sample_metadata, by = "sample_id")

#Construct covariates
covariates = dplyr::left_join(combined_expression_data$sample_metadata, peer_factors, by = "sample_id") %>%
  dplyr::mutate(sex_binary = ifelse(sex == "male", 0, 1))
rna_expressed$sample_metadata = covariates
saveRDS(rna_expressed, "results/SL1344/combined_expression_data_covariates.rds")
