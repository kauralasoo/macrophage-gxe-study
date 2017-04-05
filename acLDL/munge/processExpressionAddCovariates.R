library("devtools")
library("plyr")
library("dplyr")
library("SummarizedExperiment")
load_all("../seqUtils/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/acLDL/acLDL_combined_expression_data.rds")

#Keep only acLDL and Ctrl conditions
expression_data = extractConditionFromExpressionList(c("Ctrl", "AcLDL"), combined_expression_data)

#Filter genes by expression level
mean_expression = calculateMean(expression_data$cqn, as.data.frame(expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_expression, 1, max) > 0.5))
not_Y_genes = dplyr::filter(expression_data$gene_metadata, !(chr %in% c("MT","Y")))$gene_id
keep_genes = intersect(expressed_genes, not_Y_genes)
rna_expressed = extractGenesFromExpressionList(expression_data, keep_genes)

#Extract separate lists for each condition
condition_names = idVectorToList(c("Ctrl","AcLDL"))
rna_conditions = lapply(condition_names, extractConditionFromExpressionList, rna_expressed)

#Rename column names to genotype ids
rna_conditions_renamed = lapply(rna_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#### Set up covariates ####
#Save expression data for PEER and run PEER outside of R
cqn_list = lapply(rna_conditions_renamed, function(x){x$cqn})
savePEERData(cqn_list, "results/acLDL/PEER/input/")

#Run PEER .....
#Import PEER residuals
ctrl_residuals = t(readr::read_delim("results/acLDL/PEER/output/Ctrl/residuals.txt", delim = ",", col_names = FALSE))
rownames(ctrl_residuals) = rownames(rna_conditions$Ctrl$cqn)
colnames(ctrl_residuals) = colnames(rna_conditions$Ctrl$cqn)
acldl_residuals = t(readr::read_delim("results/acLDL/PEER/output/AcLDL/residuals.txt", delim = ",", col_names = FALSE))
rownames(acldl_residuals) = rownames(rna_conditions$AcLDL$cqn)
colnames(acldl_residuals) = colnames(rna_conditions$AcLDL$cqn)
peer_residuals = cbind(ctrl_residuals, acldl_residuals)
saveRDS(peer_residuals,"results/acLDL/PEER/output/PEER_residuals.rds")

#Import PEER factors into R
ctrl_factors = importPEERFactors("results/acLDL/PEER/output/Ctrl/factors.txt", rna_conditions_renamed$Ctrl$sample_metadata)
acldl_factors = importPEERFactors("results/acLDL/PEER/output/AcLDL/factors.txt", rna_conditions_renamed$AcLDL$sample_metadata)
peer_factors = rbind(ctrl_factors, acldl_factors) %>%
  dplyr::semi_join(combined_expression_data$sample_metadata, by = "sample_id")

#Construct covariates
covariates = dplyr::left_join(rna_expressed$sample_metadata, peer_factors, by = "sample_id") %>%
  dplyr::mutate(sex_binary = ifelse(sex == "male", 0, 1))
rna_expressed$sample_metadata = covariates
saveRDS(rna_expressed, "results/acLDL/acLDL_combined_expression_data_covariates.rds")

#Export the data in summarised-experiment format as well
#Construct a SummarizedExperiment object
coldata = rna_expressed$sample_metadata %>% as.data.frame()
rownames(coldata) = rna_expressed$sample_metadata$sample_id

rowdata = rna_expressed$gene_metadata %>% as.data.frame()
rownames(rowdata) = rna_expressed$gene_metadata$gene_id

se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = rna_expressed$counts, tpms = rna_expressed$tpm, 
                cqn = rna_expressed$cqn), 
  colData = coldata, 
  rowData = rowdata)
saveRDS(se, "results/acLDL/acLDL_combined_expression_data_covariates.se.rds")


