library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Extract separate lists for each condition
condition_names = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
rna_conditions = lapply(condition_names, extractConditionFromExpressionList, combined_expression_data)

#Rename column names to genotype ids
rna_conditions_renamed = lapply(rna_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#### RASQUAL ####
rasqual_input_folder = "results/SL1344/rasqual/input/"

#Count SNPs overlapping genes
snp_coords = readr::read_delim("genotypes/SL1344/imputed_20151005/imputed.86_samples.snp_coords.txt", 
                               delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))
exon_df = countSnpsOverlapingExons(rna_conditions_renamed$naive$gene_metadata, snp_coords, cis_window = 500000) %>% 
  dplyr::arrange(chromosome_name, range_start)
write.table(exon_df, file.path(rasqual_input_folder, "gene_snp_count_500kb.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

#Construct batches
chr11_batches = dplyr::filter(exon_df, chromosome_name == "11") %>%
  seqUtils::rasqualConstructGeneBatches(10)
write.table(chr11_batches, file.path(rasqual_input_folder, "chr11_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#Split genes into batches based on how many cis and feature SNPs they have
batches = rasqualOptimisedGeneBatches(exon_df, c(20,8,3,1))
write.table(batches, file.path(rasqual_input_folder, "gene_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#Export read counts
counts_list = lapply(rna_conditions_renamed, function(x){x$counts})
saveRasqualMatrices(counts_list, rasqual_input_folder, file_suffix = "expression")

#Extract sample-genotype map for each condition
sg_map = lapply(rna_conditions_renamed, function(x){ dplyr::select(x$sample_metadata, sample_id, genotype_id) })
saveFastqtlMatrices(sg_map, rasqual_input_folder, file_suffix = "sg_map", col_names = FALSE)

#Export GC-corrected library sizes
gc_library_size_list = lapply(counts_list, rasqualCalculateSampleOffsets, rna_conditions_renamed$naive$gene_metadata)
rasqualTools::saveRasqualMatrices(gc_library_size_list, rasqual_input_folder, file_suffix = "gc_library_size")

#Save peak names to disk
gene_names = rownames(counts_list[[1]])
write.table(gene_names, file.path(rasqual_input_folder, "gene_names.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

#Save covariates to disk
covariate_names = c("genotype_id", "PEER_factor_1", "PEER_factor_2", "sex_binary")
rasqual_cov_list = lapply(rna_conditions_renamed, function(x, covariate_names){
  rasqualMetadataToCovariates(x$sample_metadata[,covariate_names])
}, covariate_names)
saveRasqualMatrices(rasqual_cov_list, rasqual_input_folder, file_suffix = "covariates")

#Make new batches of failed peaks
naive_failed_batches = readr::read_csv("results/SL1344/rasqual/output/naive_500kb/naive_500kb.completed_ids.txt", col_names = "gene_id") %>%
  dplyr::anti_join(exon_df, ., by = "gene_id") %>%
  rasqualOptimisedGeneBatches(c(1,1,1,1), batch_prefix = "batch_5")
write.table(naive_failed_batches, "results/SL1344/rasqual/input/naive_failed_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

ifng_failed_batches = readr::read_csv("results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.completed_ids.txt", col_names = "gene_id") %>%
  dplyr::anti_join(exon_df, ., by = "gene_id") %>%
  rasqualOptimisedGeneBatches(c(1,1,1,1), batch_prefix = "batch_5")
write.table(ifng_failed_batches, "results/SL1344/rasqual/input/IFNg_failed_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

sl1344_failed_batches = readr::read_csv("results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.completed_ids.txt", col_names = "gene_id") %>%
  dplyr::anti_join(exon_df, ., by = "gene_id") %>%
  rasqualOptimisedGeneBatches(c(1,1,1,1), batch_prefix = "batch_5")
write.table(sl1344_failed_batches, "results/SL1344/rasqual/input/SL1344_failed_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

ifng_sl1344_failed_batches = readr::read_csv("results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.completed_ids.txt", col_names = "gene_id") %>%
  dplyr::anti_join(exon_df, ., by = "gene_id") %>%
  rasqualOptimisedGeneBatches(c(1,1,1,1), batch_prefix = "batch_5")
write.table(ifng_sl1344_failed_batches, "results/SL1344/rasqual/input/IFNg_SL1344_failed_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



### FastQTL ####
fastqtl_input_folder = "results/SL1344/fastqtl/input/"
#exportDataForFastQTL(rna_conditions_renamed, fastqtl_input_folder, n_chunks = 200)

#### Export data for FastQTL ####
fastqtl_genepos = constructFastQTLGenePos(rna_conditions_renamed$naive$gene_metadata)
cqn_list = lapply(rna_conditions_renamed, function(x){x$cqn})
fastqtl_cqn_list = lapply(cqn_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_cqn_list, fastqtl_input_folder, file_suffix = "expression_cqn")

#Save covariates
covariate_names = c("genotype_id", "PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6", "sex_binary")
covariate_list = lapply(rna_conditions_renamed, function(x, names){x$sample_metadata[,names]}, covariate_names)
fastqtl_covariates = lapply(covariate_list, fastqtlMetadataToCovariates)
saveFastqtlMatrices(fastqtl_covariates, fastqtl_input_folder, file_suffix = "covariates")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:250), n = 250)
write.table(chunks_matrix, file.path(fastqtl_input_folder, "all_chunk_table.txt"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")


### This is incorrect 
#Chr11 only
#Extract chr11 only
chr11_genes = dplyr::filter(atac_list$gene_metadata, chr == 11)$gene_id
chr11_list = lapply(atac_conditions_renamed, extractGenesFromExpressionList, chr11_genes)
fastqtl_genepos = constructFastQTLGenePos(chr11_list$naive$gene_metadata)

#Save expression data
fastqtl_cqn_list = lapply(chr11_list, function(x){x$cqn}) %>%
  lapply(., prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_cqn_list, "results/ATAC/fastqtl/input/", file_suffix = "expression")
fastqtl_tpm_list = lapply(chr11_list, function(x){x$tpm}) %>%
  lapply(., prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_tpm_list, "results/ATAC/fastqtl/input/", file_suffix = "expression_tpm")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:25), n = 25)
write.table(chunks_matrix, "results/ATAC/fastqtl/input/chunk_table.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")

#### eigenMT ####
#Export genotype data
chromosome_list = scan("macrophage-gxe-study/data/sample_lists/chromosome_list.txt", what = "char")
eigenMTExportGenotypesByChr(chromosome_list, "genotypes/SL1344/imputed_20151005/chromosomes_INFO_07/",
                            "results/SL1344/eigenMT/input/", "chr_")
#Export gene metadata
eigenMTExportGeneMetadata(combined_expression_data$gene_metadata, "results/SL1344/eigenMT/input/")
