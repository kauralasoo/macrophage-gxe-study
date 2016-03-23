library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")
library("readr")
library("rasqualTools")
load_all("~/software/rasqual/rasqualTools/")
#Import atac data list
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Extract separate lists for each condition
naive_list = extractConditionFromExpressionList("naive", atac_list)
IFNg_list = extractConditionFromExpressionList("IFNg", atac_list)
SL1344_list = extractConditionFromExpressionList("SL1344", atac_list)
IFNg_SL1344_list = extractConditionFromExpressionList("IFNg_SL1344", atac_list)
atac_conditions = list(naive = naive_list, IFNg = IFNg_list, SL1344 = SL1344_list, IFNg_SL1344 = IFNg_SL1344_list)

#Save expression data for PEER and run PEER outside of R
peer_cqn_list = lapply(atac_conditions, function(x){x$cqn})
savePEERData(peer_cqn_list, "results/ATAC/PEER/input/")

peer_tpm_list = lapply(atac_conditions, function(x){log(x$tpm + 0.01, 2)})
savePEERData(peer_cqn_list, "results/ATAC/PEER/input/",file_suffix = "exprs_tpm")

#Import PEER factors back into R
naive_peer = importPEERFactors("results/ATAC/PEER/output/naive_10/factors.txt", atac_conditions$naive$sample_metadata)
IFNg_peer = importPEERFactors("results/ATAC/PEER/output/IFNg_10/factors.txt", atac_conditions$IFNg$sample_metadata)
SL1344_peer = importPEERFactors("results/ATAC/PEER/output/SL1344_10/factors.txt", atac_conditions$SL1344$sample_metadata)
IFNg_SL1344_peer = importPEERFactors("results/ATAC/PEER/output/IFNg_SL1344_10/factors.txt", atac_conditions$IFNg_SL1344$sample_metadata)
peer_factors = rbind(naive_peer, IFNg_peer, SL1344_peer, IFNg_SL1344_peer)
sample_ids = peer_factors$sample_id
covariates = t(peer_factors[,-1])
colnames(covariates)= sample_ids
covariates = covariates[,colnames(head(atac_list$counts))]

#Use PCA to construct covariates for QTL mapping
naive_PCA = performPCA(log(atac_conditions$naive$tpm + 0.01, 2), atac_conditions$naive$sample_metadata, n_pcs = 5)
IFNg_PCA = performPCA(log(atac_conditions$IFNg$tpm + 0.01, 2), atac_conditions$IFNg$sample_metadata, n_pcs = 5)
SL1344_PCA = performPCA(log(atac_conditions$SL1344$tpm + 0.01, 2), atac_conditions$SL1344$sample_metadata, n_pcs = 5)
IFNg_SL1344_PCA = performPCA(log(atac_conditions$IFNg_SL1344$tpm + 0.01, 2), atac_conditions$IFNg_SL1344$sample_metadata, n_pcs = 5)

#Plot variance explained
plot(naive_PCA$var_exp)
plot(IFNg_PCA$var_exp)
plot(SL1344_PCA$var_exp)
plot(IFNg_SL1344_PCA$var_exp)

#Add PCs into the sample metadata df
covariates = rbind(naive_PCA$pca_matrix, IFNg_PCA$pca_matrix, SL1344_PCA$pca_matrix, IFNg_SL1344_PCA$pca_matrix) %>% 
  tbl_df() %>%
  dplyr::mutate(sex_binary = ifelse(sex == "male",0,1))
new_sample_metadata = dplyr::select(atac_list$sample_metadata, sample_id) %>% dplyr::left_join(covariates, by = "sample_id")
atac_list$sample_metadata = new_sample_metadata

#Extract separate lists for each condition
naive_list = extractConditionFromExpressionList("naive", atac_list)
IFNg_list = extractConditionFromExpressionList("IFNg", atac_list)
SL1344_list = extractConditionFromExpressionList("SL1344", atac_list)
IFNg_SL1344_list = extractConditionFromExpressionList("IFNg_SL1344", atac_list)
atac_conditions = list(naive = naive_list, IFNg = IFNg_list, SL1344 = SL1344_list, IFNg_SL1344 = IFNg_SL1344_list)

#Rename column names to genotype ids
atac_conditions_renamed = lapply(atac_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#### Export data for FastQTL ####
tpm_list = lapply(atac_conditions_renamed, function(x){x$tpm})
fastqtl_genepos = constructFastQTLGenePos(naive_list$gene_metadata)
fastql_tpm_list = lapply(tpm_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastql_tpm_list, "results/ATAC/fastqtl/input/", file_suffix = "tpm_exp_full")

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

#Save covariates
fastqtl_cov_list = lapply(chr11_list, function(x){x$covariates}) %>%
  lapply(., prepareFastqtlCovariates, 4:5)
saveFastqtlMatrices(fastqtl_cov_list, "results/ATAC/fastqtl/input/", file_suffix = "covariates_pc2")

fastqtl_cov_list = lapply(chr11_list, function(x){x$covariates}) %>%
  lapply(., prepareFastqtlCovariates, 1:6)
saveFastqtlMatrices(fastqtl_cov_list, "results/ATAC/fastqtl/input/", file_suffix = "covariates_pc3")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:25), n = 25)
write.table(chunks_matrix, "results/ATAC/fastqtl/input/chunk_table.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")

chunks_matrix = data.frame(chunk = seq(1:200), n = 200)
write.table(chunks_matrix, "results/ATAC/fastqtl/input/all_chunk_table.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")


#### export data for RASQUAL ####
counts_list = lapply(atac_conditions_renamed, function(x){x$counts})
saveRasqualMatrices(counts_list, "results/ATAC/rasqual/input/", file_suffix = "expression")

#Extract sample-genotype map for each condition
sg_map = lapply(atac_conditions_renamed, function(x){ dplyr::select(x$sample_metadata, sample_id, genotype_id) })
saveFastqtlMatrices(sg_map,"results/ATAC/rasqual/input/", file_suffix = "sg_map", col_names = FALSE)

#Count overlaps between SNPs and peaks
snp_coords = readr::read_delim("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.snp_coords.txt", 
                               delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))
peak_snp_count_100kb = rasqualTools::countSnpsOverlapingPeaks(atac_list$gene_metadata, snp_coords, cis_window = 100000)
write.table(peak_snp_count_100kb, "results/ATAC/rasqual/input/peak_snp_count_100kb.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Export GC-corrected library sizes
gc_library_size_list = lapply(counts_list, rasqualCalculateSampleOffsets, atac_conditions_renamed[[1]]$gene_metadata)
rasqualTools::saveRasqualMatrices(gc_library_size_list, "results/ATAC/rasqual/input/", file_suffix = "gc_library_size")

#Save peak names to disk
peak_names = rownames(atac_list$counts)
write.table(peak_names, "results/ATAC/rasqual/input/peak_names.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Save covariates to disk
covariate_names = c("genotype_id", "PC1", "PC2", "PC3", "sex_binary")
rasqual_cov_list = lapply(atac_conditions_renamed, function(x, covariate_names){
  rasqualMetadataToCovariates(x$sample_metadata[,covariate_names])
  }, covariate_names)
saveRasqualMatrices(rasqual_cov_list, "results/ATAC/rasqual/input/", file_suffix = "covariates")

#Construct batches
chr11_batches = dplyr::filter(atac_list$gene_metadata, chr == 11) %>%
  seqUtils::rasqualConstructGeneBatches(70)
write.table(chr11_batches, "results/ATAC/rasqual/input/peak_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#Construct batches for all genes
atac_batches = rasqualTools::rasqualOptimisedGeneBatches(peak_snp_count_100kb, batch_sizes = c(150, 20, 3,1))
write.table(atac_batches, "results/ATAC/rasqual/input/all_peak_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#Make new batches of failed peaks
naive_failed_batches = readr::read_csv("results/ATAC/rasqual/output/naive_100kb/naive_100kb.completed_ids.txt", col_names = "gene_id") %>%
  dplyr::anti_join(peak_snp_count_100kb, ., by = "gene_id") %>%
  rasqualOptimisedGeneBatches(c(10,10,1,1))
write.table(naive_failed_batches, "results/ATAC/rasqual/input/naive_failed_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

IFNg_failed_batches = readr::read_csv("results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.completed_ids.txt", col_names = "gene_id") %>%
  dplyr::anti_join(peak_snp_count_100kb, ., by = "gene_id") %>%
  rasqualOptimisedGeneBatches(c(10,10,1,1))
write.table(IFNg_failed_batches, "results/ATAC/rasqual/input/IFNg_failed_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

SL1344_failed_batches = readr::read_csv("results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.completed_ids.txt", col_names = "gene_id") %>%
  dplyr::anti_join(peak_snp_count_100kb, ., by = "gene_id") %>%
  rasqualOptimisedGeneBatches(c(10,10,1,1))
write.table(SL1344_failed_batches, "results/ATAC/rasqual/input/SL1344_failed_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

IFNg_SL1344_failed_batches = readr::read_csv("results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.completed_ids.txt", col_names = "gene_id") %>%
  dplyr::anti_join(peak_snp_count_100kb, ., by = "gene_id") %>%
  rasqualOptimisedGeneBatches(c(10,10,1,1), batch_prefix = "batch_5")
write.table(IFNg_SL1344_failed_batches, "results/ATAC/rasqual/input/IFNg_SL1344_failed_batches.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#Export ATAC peak positions for eigenMT
eigenMTExportGeneMetadata(atac_list$gene_metadata, "results/ATAC/eigenMT/input")
