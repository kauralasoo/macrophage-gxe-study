library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")

#Import atac data list
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata = dplyr::mutate(atac_list$sample_metadata, genotype_id = ifelse(genotype_id == "HPSI1213i-nusw_2", "HPSI1213i-nusw_1", genotype_id)) #Fix this temportary bug in genotypes vcf file

#Extract separate lists for each condition
naive_list = extractConditionFromExpressionList(atac_list, "naive")
IFNg_list = extractConditionFromExpressionList(atac_list, "IFNg")
SL1344_list = extractConditionFromExpressionList(atac_list, "SL1344")
IFNg_SL1344_list = extractConditionFromExpressionList(atac_list, "IFNg_SL1344")
atac_conditions = list(naive = naive_list, IFNg = IFNg_list, SL1344 = SL1344_list, IFNg_SL1344 = IFNg_SL1344_list)

#Rename column names to genotype ids
atac_conditions_renamed = lapply(atac_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#Save expression data for PEER and run PEER outside of R
peer_cqn_list = lapply(atac_conditions_renamed, function(x){x$cqn})
savePEERData(peer_cqn_list, "results/ATAC/PEER/input/")

peer_tpm_list = lapply(atac_conditions_renamed, function(x){x$tpm})
savePEERData(peer_cqn_list, "results/ATAC/PEER/input/",file_suffix = "exprs_tpm")

#### Export data for FastQTL ####
cqn_list = lapply(atac_conditions_renamed, function(x){x$cqn})
fastqtl_genepos = constructFastQTLGenePos(naive_list$gene_metadata)
fastql_cqn_list = lapply(cqn_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastql_cqn_list, "results/ATAC/fastqtl/input/", file_suffix = "expression")

#Extract chr21 only
chr21_genes = dplyr::filter(atac_list$gene_metadata, chr == 21)$gene_id
chr21_list = lapply(atac_conditions_renamed, extractGenesFromExpressionList, chr21_genes)
fastqtl_genepos = constructFastQTLGenePos(chr21_list$naive$gene_metadata)
fastqtl_cqn_list = lapply(chr21_list, function(x){x$cqn}) %>%
  lapply(., prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_cqn_list, "results/ATAC/fastqtl/input/", file_suffix = "expression")
fastqtl_tpm_list = lapply(chr21_list, function(x){x$tpm}) %>%
  lapply(., prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_tpm_list, "results/ATAC/fastqtl/input/", file_suffix = "expression_tpm")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:25), n = 25)
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
peak_snp_count_50kb = seqUtils::countSnpsOverlapingPeaks(atac_list$gene_metadata, snp_coords, cis_window = 50000)
write.table(peak_snp_count, "results/ATAC/rasqual/input/peak_snp_count_500kb.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(peak_snp_count_2kb, "results/ATAC/rasqual/input/peak_snp_count_2kb.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(peak_snp_count_50kb, "results/ATAC/rasqual/input/peak_snp_count_50kb.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Extract sample offsests
norm_factors_list = lapply(atac_conditions_renamed, rasqualSizeFactorsMatrix, "norm_factor")
saveRasqualMatrices(norm_factors_list, "results/ATAC/rasqual/input/", file_suffix = "offsets")
library_size_list = lapply(atac_conditions_renamed, rasqualSizeFactorsMatrix, "library_size")
saveRasqualMatrices(library_size_list, "results/ATAC/rasqual/input/", file_suffix = "library_size")

#Save peak names to disk
peak_names = rownames(atac_list$counts)
write.table(peak_names, "results/ATAC/rasqual/input/peak_names.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Construct batches
chr21_batches = dplyr::filter(atac_list$gene_metadata, chr == 21) %>% 
  dplyr::select(gene_id) %>%
  dplyr::mutate(batch_number = splitIntoBatches(length(gene_id), 35)) %>%
  dplyr::mutate(batch_id = paste("batch", batch_number, sep = "_")) %>% 
  dplyr::group_by(batch_id) %>%
  dplyr::summarize(gene_ids = paste(gene_id, collapse = ","))
write.table(chr21_batches, "results/ATAC/rasqual/input/peak_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

