library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")

#Import atac data list
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Extract separate lists for each condition
naive_list = extractConditionFromExpressionList(atac_list, "naive")
IFNg_list = extractConditionFromExpressionList(atac_list, "IFNg")
SL1344_list = extractConditionFromExpressionList(atac_list, "SL1344")
IFNg_SL1344_list = extractConditionFromExpressionList(atac_list, "IFNg_SL1344")
atac_conditions = list(naive = naive_list, IFNg = IFNg_list, SL1344 = SL1344_list, IFNg_SL1344 = IFNg_SL1344_list)

#Rename column names to genotype ids
atac_conditions_renamed = lapply(atac_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#### Export data for FastQTL ####
cqn_list = lapply(atac_conditions_renamed, function(x){x$cqn})
fastqtl_genepos = constructFastQTLGenePos(naive_list$gene_metadata)
fastql_cqn_list = lapply(cqn_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastql_cqn_list, "results/ATAC/fastqtl/input/", file_suffix = "expression")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:200), n = 200)
write.table(chunks_matrix, "results/ATAC/fastqtl/input/chunk_table.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")


#### export data for RASQUAL ####
counts_list = lapply(atac_conditions_renamed, function(x){x$counts})
saveRasqualMatrices(counts_list, "results/ATAC/rasqual/input/", file_suffix = "expression")

#Extract sample-genotype map
sg_map = lapply(atac_conditions_renamed, function(x){ dplyr::select(x$sample_metadata, sample_id, genotype_id) })
saveFastqtlMatrices(sg_map,"results/ATAC/rasqual/input/", file_suffix = "sg_map", col_names = FALSE)