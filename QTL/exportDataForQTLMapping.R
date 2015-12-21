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

#Extract sample-genotype map for each condition
sg_map = lapply(atac_conditions_renamed, function(x){ dplyr::select(x$sample_metadata, sample_id, genotype_id) })
saveFastqtlMatrices(sg_map,"results/ATAC/rasqual/input/", file_suffix = "sg_map", col_names = FALSE)

#Count overlaps between SNPs and peaks
snp_coords = read.table("results/ATAC/rasqual/input/imputed.69_samples.snp_coords.txt", stringsAsFactors = FALSE)
colnames(snp_coords) = c("chr", "pos", "snp_id")
peak_snp_count = seqUtils::countSnpsOverlapingPeaks(atac_list$gene_metadata, snp_coords, cis_window = 500000)
peak_snp_count_2kb = seqUtils::countSnpsOverlapingPeaks(atac_list$gene_metadata, snp_coords, cis_window = 2000)
write.table(peak_snp_count, "results/ATAC/rasqual/input/peak_snp_count_500kb.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(peak_snp_count_2kb, "results/ATAC/rasqual/input/peak_snp_count_2kb.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Extract sample offsests
norm_factors_list = lapply(atac_conditions_renamed, rasqualSizeFactorsMatrix, "norm_factor")
saveRasqualMatrices(norm_factors_list, "results/ATAC/rasqual/input/", file_suffix = "offsets")

#Save peak names to disk
peak_names = rownames(atac_list$counts)
write.table(peak_names, "results/ATAC/rasqual/input/peak_names.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

