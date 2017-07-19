library("dplyr")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")


#Export full eQTL dataset
#Import SummarizedExperiment
se = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")
event_metadata = rowData(se)
sample_metadata = colData(se) %>% tbl_df2()
unique_genes = unique(event_metadata$probe_id)

#Remove events on X and Y chromosomes
event_dataset = se[event_metadata[!event_metadata$chr %in% c("X","Y","MT"),]$probe_id,]

#Extract lists for each condition
condition_list = idVectorToList(c("CD14","IFN","LPS2","LPS24"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,event_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Expression level filtering (NA)
expressed_dataset = event_dataset

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(expressed_dataset)) %>% 
  dplyr::filter(probe_id %in% rownames(expressed_dataset)) %>%
  dplyr::mutate(transcript_id = probe_id, start = gene_start, end = gene_end) %>%
  constructQTLtoolsGenePos()
output_path = "processed/Fairfax/qtltools/input/full/"

#Extract cqn matrices
exprs_list = purrr::map(event_conditions_renamed, ~assays(.)$exprs)
qtltools_cqn_list = purrr::map(exprs_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_cqn_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = purrr::map(event_conditions_renamed, ~tbl_df2(colData(.)))
covariates_list = purrr::map2(exprs_list, sample_meta, 
                             ~performPCA(.x, .y, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")


#Export dataset of samples shared between all conditions (shared)
#Import SummarizedExperiment
se = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")
event_metadata = rowData(se)
sample_metadata = colData(se) %>% tbl_df2()
unique_genes = unique(event_metadata$probe_id)

#Remove events on X and Y chromosomes
event_dataset = se[event_metadata[!event_metadata$chr %in% c("X","Y","MT"),]$probe_id,]

#Keep only those donors that are present in all conditions
event_dataset = event_dataset[,dplyr::filter(sample_metadata, present_in_all)$sample_id]

#Extract lists for each condition
condition_list = idVectorToList(c("CD14","IFN","LPS2","LPS24"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,event_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Expression level filtering (NA)
expressed_dataset = event_dataset

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(expressed_dataset)) %>% 
  dplyr::filter(probe_id %in% rownames(expressed_dataset)) %>%
  dplyr::mutate(transcript_id = probe_id, start = gene_start, end = gene_end) %>%
  constructQTLtoolsGenePos()
output_path = "processed/Fairfax/qtltools/input/shared/"

#Extract cqn matrices
exprs_list = purrr::map(event_conditions_renamed, ~assays(.)$exprs)
qtltools_cqn_list = purrr::map(exprs_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_cqn_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = purrr::map(event_conditions_renamed, ~tbl_df2(colData(.)))
covariates_list = purrr::map2(exprs_list, sample_meta, 
                              ~performPCA(.x, .y, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                                dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                                fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")





#Export dataset for a random subset of individuals (84)
#Import SummarizedExperiment
se = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")
event_metadata = rowData(se)
sample_metadata = colData(se) %>% tbl_df2()
unique_genes = unique(event_metadata$probe_id)

#Remove events on X and Y chromosomes
event_dataset = se[event_metadata[!event_metadata$chr %in% c("X","Y","MT"),]$probe_id,]

#Keep only those donors that are present in all conditions and choose a random sample
set.seed(42)
shared_samples = dplyr::filter(sample_metadata, present_in_all)
shared_donors = unique(shared_samples$donor_id)
selected_samples = dplyr::filter(shared_samples, donor_id %in% sample(shared_donors, 84))
event_dataset = event_dataset[,selected_samples$sample_id]

#Extract lists for each condition
condition_list = idVectorToList(c("CD14","IFN","LPS2","LPS24"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,event_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Expression level filtering (NA)
expressed_dataset = event_dataset

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(expressed_dataset)) %>% 
  dplyr::filter(probe_id %in% rownames(expressed_dataset)) %>%
  dplyr::mutate(transcript_id = probe_id, start = gene_start, end = gene_end) %>%
  constructQTLtoolsGenePos()
output_path = "processed/Fairfax/qtltools/input/shared_84/"

#Extract cqn matrices
exprs_list = purrr::map(event_conditions_renamed, ~assays(.)$exprs)
qtltools_cqn_list = purrr::map(exprs_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_cqn_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = purrr::map(event_conditions_renamed, ~tbl_df2(colData(.)))
covariates_list = purrr::map2(exprs_list, sample_meta, 
                              ~performPCA(.x, .y, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                                dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                                fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")


