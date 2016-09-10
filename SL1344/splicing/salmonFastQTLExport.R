library("DRIMSeq")
library("SummarizedExperiment")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import event dataset
se_ensembl = readRDS("results/SL1344/combined_reviseAnnotations_transcript_quants.rds")
event_metadata = rowData(se_ensembl)
unique_genes = unique(event_metadata$ensembl_gene_id)

#Calculate mean TPM in each condition
gene_transcript_map = tbl_df2(rowData(se_ensembl)) %>% dplyr::select(transcript_id, gene_id)
tpm_matrix = assays(se_ensembl)$tpm
tpm_df = dplyr::mutate(as.data.frame(tpm_matrix), transcript_id = rownames(tpm_matrix)) %>%
  tbl_df() %>%
  dplyr::left_join(gene_transcript_map, by = "transcript_id") %>%
  dplyr::select(-transcript_id) %>%
  dplyr::select(gene_id, everything())

#Calculate total expression per gene
mean_gene_expression = purrr::slice_rows(tpm_df, "gene_id") %>% 
  purrr::by_slice(~colSums(.) %>% 
                    t() %>% 
                    as.data.frame(), .collate = "rows")

#Make matrix of gene expression values for each transcript
tx_gene_expression = dplyr::left_join(gene_transcript_map, mean_gene_expression, by = "gene_id")
tx_gene_matrix = dplyr::select(tx_gene_expression, -gene_id, -transcript_id) %>% as.matrix()
rownames(tx_gene_matrix) = tx_gene_expression$transcript_id

#calculate TPM ratios and add them to the SummarizedExperiment
tpm_ratios = tpm_matrix/tx_gene_matrix[rownames(tpm_matrix),]
assays(se_ensembl)$tpm_ratios = tpm_ratios

#Mean transript expression
mean_tx_expression = calculateMean(tpm_matrix, colData(se_ensembl), "condition_name")
expressed_transcripts = names(which(apply(mean_tx_expression, 1, max) > 1))

#Mean gene expression
gene_expression_matrix = tibbleToNamedMatrix(mean_gene_expression, row_names = "gene_id")
mean_gene_exp = calculateMean(gene_expression_matrix, colData(se_ensembl), "condition_name")
expressed_genes = names(which(apply(mean_gene_exp, 1, max) > 1))

#Filter transcript ratios
tpm_ratios_filtered = se_ensembl[dplyr::filter(gene_transcript_map, gene_id %in% expressed_genes)$transcript_id,]

#Extract lists for each condition
condition_list = idVectorToList(c("naive","IFNg", "SL1344","IFNg_SL1344"))
conditions_se = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,tpm_ratios_filtered))

#Rename columns
renamed_se_list = purrr::map(conditions_se, function(x){
  colnames(x) = x$genotype_id
  return(x)
  })

#### Export data for FastQTL ####
fastqtl_genepos = tbl_df2(rowData(se_ensembl)) %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::mutate(start = min(transcript_start), end = max(transcript_end)) %>% 
  dplyr::ungroup() %>% 
  dplyr::transmute(gene_id = transcript_id, start, end, chr) %>%
  dplyr::filter(gene_id %in% rownames(tpm_ratios_filtered)) %>%
  constructFastQTLGenePos()

#Extract matrices
prop_list = purrr::map(renamed_se_list, ~assays(.)$tpm_ratios)

#Quantile normalise and remove NAs
normalised_list = purrr::map(prop_list, ~replaceNAsWithRowMeans(.)) %>%
  purrr::map(~t(quantileNormaliseMatrix(t(.))))
fastqtl_norm_prop_list = purrr::map(normalised_list, prepareFastqtlMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_norm_prop_list, "results/SL1344/salmon/fastqtl_input/", file_suffix = "norm_prop")

#Calculate principal components
sample_meta = tbl_df2(colData(renamed_se_list$naive)) %>% 
  dplyr::rename(sample_id2 = sample_id) %>%
  dplyr::rename(sample_id = genotype_id)

#Perform PCA on each dataset and construct covariates
pca_list = purrr::map(normalised_list, ~performPCA(., sample_meta, n_pcs = 5)$pca_matrix %>% 
                        dplyr::rename(genotype_id = sample_id) %>%
                        dplyr::rename(sample_id = sample_id2))
selected_covariates = purrr::map(pca_list, ~dplyr::select(.,genotype_id, PC1, PC2, PC3, PC4, PC5) %>%
                                   fastqtlMetadataToCovariates())
saveFastqtlMatrices(selected_covariates, "results/SL1344/salmon/fastqtl_input/", file_suffix = "covariates_prop")

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:250), n = 250)
write.table(chunks_matrix, "results/SL1344/salmon/fastqtl_input/all_chunk_table.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")




