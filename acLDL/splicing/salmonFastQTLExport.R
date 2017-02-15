library("SummarizedExperiment")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import event dataset
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
event_metadata = rowData(se_ensembl)
unique_genes = unique(event_metadata$ensembl_gene_id)

#Remove events on X and Y chromosomes
event_dataset = se_ensembl[event_metadata[!event_metadata$chr %in% c("X","Y"),]$transcript_id,]

#Extract lists for each condition
condition_list = idVectorToList(c("Ctrl","AcLDL"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,event_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(event_dataset)) %>% 
  dplyr::transmute(gene_id = transcript_id, start, end, chr) %>%
  dplyr::filter(gene_id %in% rownames(event_dataset)) %>%
  constructFastQTLGenePos()

#Extract ratio matrices
prop_list = purrr::map(event_conditions_renamed, ~assays(.)$tpm_ratios)

#Quantile normalise and remove NAs
normalised_list = purrr::map(prop_list, ~replaceNAsWithRowMeans(.) %>% quantileNormaliseRows())
fastqtl_norm_prop_list = purrr::map(normalised_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_norm_prop_list, "results/acLDL/fastqtl_splicing/ensembl_87/", file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$Ctrl))
covariates_list = purrr::map(normalised_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, "results/acLDL/fastqtl_splicing/ensembl_87/", file_suffix = "covariates_prop")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:250), n = 250)
write.table(chunks_matrix, "results/acLDL/fastqtl_splicing/ensembl_87/all_chunk_table.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")

