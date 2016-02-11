library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")

#Import data
combined_expression_data = readRDS("results/acLDL/acLDL_combined_expression_data_covariates.rds")

#Extract separate lists for each condition
condition_names = idVectorToList(c("Ctrl","AcLDL"))
rna_conditions = lapply(condition_names, extractConditionFromExpressionList, combined_expression_data)

#Rename column names to genotype ids
rna_conditions_renamed = lapply(rna_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#### RASQUAL ####
rasqual_input_folder = "results/acLDL/rasqual/input/"
exportDataForRasqual(rna_conditions_renamed, rasqual_input_folder)

#Count SNPs overlapping genes
snp_coords = readr::read_delim("genotypes/acLDL/imputed_20151005/imputed.70_samples.snp_coords.txt", 
                               delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))
exon_df = countSnpsOverlapingExons(rna_conditions_renamed$Ctrl$gene_metadata, snp_coords, cis_window = 500000) %>% 
  dplyr::arrange(chromosome_name, range_start)
write.table(exon_df, file.path(rasqual_input_folder, "gene_snp_count_500kb.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

exon_df = countSnpsOverlapingExons(rna_conditions_renamed$Ctrl$gene_metadata, snp_coords, cis_window = 100000) %>% 
  dplyr::arrange(chromosome_name, range_start)
write.table(exon_df, file.path(rasqual_input_folder, "gene_snp_count_100kb.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

#Construct batches
chr11_batches = dplyr::filter(exon_df, chromosome_name == "11") %>%
  seqUtils::rasqualConstructGeneBatches(10)
write.table(chr11_batches, file.path(rasqual_input_folder, "chr11_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#Split genes into batches based on how many cis and feature SNPs they have
batches = rasqualOptimisedGeneBatches(exon_df, c(20,8,3,1))
write.table(batches, file.path(rasqual_input_folder, "gene_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


### FastQTL ####
#Export expression data for fastQTL analysis
eqtl_dataset = readRDS("results/acLDL/acLDL_eqtl_data_list.rds")

#Extract data from list
donor_genotype_map = dplyr::filter(eqtl_dataset$sample_metadata, condition == "Ctrl") %>% dplyr::select(donor, genotype_id)
genepos = dplyr::select(eqtl_dataset$genepos, chr, left, right, geneid)
colnames(genepos)[1] = "#chr"

#Make an updated list of covariates
fastqtl_covariates_list = lapply(eqtl_dataset$covariates_list, function(cov_mat, donor_genotype_map){
  res = extractSubset(donor_genotype_map, as.data.frame(cov_mat), "donor", "genotype_id") %>% 
    dplyr::mutate(id = rownames(cov_mat)) %>% 
    dplyr::select(id, everything())
  return(res)
}, donor_genotype_map) %>%
  lapply(., function(x){x[1:7,]}) #Keep the first 7 covariates

#Make an updated list of gene expression
fastqtl_expression_list = lapply(eqtl_dataset$exprs_cqn_list, function(exp_mat, donor_genotype_map, genepos){
  res = extractSubset(donor_genotype_map, as.data.frame(exp_mat), "donor","genotype_id") %>%
    dplyr::mutate(geneid = rownames(exp_mat)) %>%
    dplyr::select(geneid, everything()) %>%
    dplyr::left_join(genepos, ., by = "geneid")
  return(res)
}, donor_genotype_map, genepos)

#Save matrices to disk
saveFastqtlMatrices(fastqtl_covariates_list, "results/acLDL/fastqtl/input/", file_suffix = "covariates")
saveFastqtlMatrices(fastqtl_expression_list, "results/acLDL/fastqtl/input/", file_suffix = "expression")
write.table(donor_genotype_map$genotype_id, "results/acLDL/fastqtl/input/genotype_list.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:200), n = 200)
write.table(chunks_matrix, "results/acLDL/fastqtl/input/chunk_table.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")



