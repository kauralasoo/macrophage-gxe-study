library("SummarizedExperiment")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import event dataset
se_ensembl = readRDS("results/SL1344/combined_reviseAnnotations_transcript_quants.rds")
event_metadata = rowData(se_ensembl)
unique_genes = unique(event_metadata$ensembl_gene_id)

#Remove events on X and Y chromosomes
event_dataset = se_ensembl[event_metadata[!event_metadata$chr %in% c("X","Y"),]$transcript_id,]

#Extract lists for each condition
condition_list = idVectorToList(c("naive","IFNg", "SL1344","IFNg_SL1344"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,event_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(event_dataset)) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(start = min(transcript_start), end = max(transcript_end)) %>% 
  dplyr::ungroup() %>% 
  dplyr::transmute(gene_id = transcript_id, start, end, chr) %>%
  dplyr::filter(gene_id %in% rownames(event_dataset)) %>%
  constructFastQTLGenePos()

#Extract ratio matrices
prop_list = purrr::map(event_conditions_renamed, ~assays(.)$tpm_ratios)

#Quantile normalise and remove NAs
normalised_list = purrr::map(prop_list, ~replaceNAsWithRowMeans(.) %>% quantileNormaliseRows())
fastqtl_norm_prop_list = purrr::map(normalised_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_norm_prop_list, "results/SL1344/salmon/fastqtl_input/", file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, "results/SL1344/salmon/fastqtl_input/", file_suffix = "covariates_prop")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:250), n = 250)
write.table(chunks_matrix, "results/SL1344/salmon/fastqtl_input/all_chunk_table.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")




#Export complete Ensmebl 85 transcript ratios for FastQTL
se_ensembl = readRDS("results/SL1344/combined_ensembl85_transcript_quants.rds")
event_metadata = tbl_df2(rowData(se_ensembl))
multi_transcript_genes = names(which(table(event_metadata$gene_id) > 1))

#Remove single transcript genes and other ivalid transcripts
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9")
valid_gene_biotypes = c("lincRNA","protein_coding")
valid_transcript_biotypes = c("lincRNA", "protein_coding")
filtered_transcripts = dplyr::filter(event_metadata, gene_id %in% multi_transcript_genes, 
              gene_biotype %in% valid_gene_biotypes, 
              chromosome_name %in% valid_chromosomes, 
              transcript_biotype %in% valid_transcript_biotypes)
filtered_events = se_ensembl[filtered_transcripts$transcript_id,]


#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(filtered_events)) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(start = min(transcript_start), end = max(transcript_end)) %>% 
  dplyr::ungroup() %>% 
  dplyr::transmute(gene_id = transcript_id, start, end, chr = chromosome_name) %>%
  dplyr::filter(gene_id %in% rownames(filtered_events)) %>%
  constructFastQTLGenePos()

#Extract lists for each condition
condition_list = idVectorToList(c("naive","IFNg", "SL1344","IFNg_SL1344"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,filtered_events))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Extract ratio matrices
prop_list = purrr::map(event_conditions_renamed, ~assays(.)$tpm_ratios)

#Quantile normalise and remove NAs
normalised_list = purrr::map(prop_list, ~replaceNAsWithRowMeans(.)) %>%
  purrr::map(~t(quantileNormaliseMatrix(t(.))))
fastqtl_norm_prop_list = purrr::map(normalised_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_norm_prop_list, "results/SL1344/salmon_ensembl85/fastqtl_input/", file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, "results/SL1344/salmon_ensembl85/fastqtl_input/", file_suffix = "covariates_prop")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:250), n = 250)
write.table(chunks_matrix, "results/SL1344/salmon_ensembl85/fastqtl_input/all_chunk_table.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")


#Export complete Ensmebl 85 transcript ratios for FastQTL without any filtering (except chromosomes)
se_ensembl = readRDS("results/SL1344/combined_ensembl85_transcript_quants.rds")
event_metadata = tbl_df2(rowData(se_ensembl))
multi_transcript_genes = names(which(table(event_metadata$gene_id) > 1))

#Remove single transcript genes and other ivalid transcripts
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9")
filtered_transcripts = dplyr::filter(event_metadata, gene_id %in% multi_transcript_genes, 
                                     chromosome_name %in% valid_chromosomes)
filtered_events = se_ensembl[filtered_transcripts$transcript_id,]


#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(filtered_events)) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(start = min(transcript_start), end = max(transcript_end)) %>% 
  dplyr::ungroup() %>% 
  dplyr::transmute(gene_id = transcript_id, start, end, chr = chromosome_name) %>%
  dplyr::filter(gene_id %in% rownames(filtered_events)) %>%
  constructFastQTLGenePos()

#Extract lists for each condition
condition_list = idVectorToList(c("naive","IFNg", "SL1344","IFNg_SL1344"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,filtered_events))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Extract ratio matrices
prop_list = purrr::map(event_conditions_renamed, ~assays(.)$tpm_ratios)

#Quantile normalise and remove NAs
normalised_list = purrr::map(prop_list, ~replaceNAsWithRowMeans(.)) %>%
  purrr::map(~t(quantileNormaliseMatrix(t(.))))
fastqtl_norm_prop_list = purrr::map(normalised_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_norm_prop_list, "results/SL1344/salmon_ensembl85_full/fastqtl_input/", file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, "results/SL1344/salmon_ensembl85_full/fastqtl_input/", file_suffix = "covariates_prop")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:400), n = 400)
write.table(chunks_matrix, "results/SL1344/salmon_ensembl85_full/fastqtl_input/all_chunk_table.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")



