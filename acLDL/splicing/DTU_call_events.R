library("dplyr")
library("DRIMSeq")
library("devtools")
library("optparse")
library("SummarizedExperiment")
load_all("../seqUtils/")

extractDTUResults <- function(result_list, se){
  
  #Extract pvalues
  pvalues_df = purrr::map_df(result_list, ~DRIMSeq::results(.)) %>% 
    dplyr::mutate(adj_pvalue = p.adjust(pvalue, method = "fdr")) %>%
    dplyr::rename(p_nominal = pvalue, p_fdr = adj_pvalue) %>%
    dplyr::tbl_df()
  
  #Calculate mean proportions
  proportions_df = purrr::map_df(result_list, ~DRIMSeq::proportions(.))
  mean_proportions = calculateMeanProportions(proportions_df, se)
  
  #Calculate max diff within cluster
  max_diff = dplyr::mutate(mean_proportions, Diff = AcLDL - Ctrl) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(max_diff = max(abs(Diff)))
  
  results = dplyr::left_join(pvalues_df, max_diff, by = "gene_id")
  return(results)
}

#Import SummarizedExperiments for all phenotypes
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
se_revised = readRDS("results/acLDL/acLDL_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/acLDL/acLDL_leafcutter_counts.rds")
se_list = list(ensembl_87 = se_ensembl, revisedAnnotation = se_revised, leafcutter = se_leafcutter)

#Extract sample metadata
sample_meta_list = purrr::map(se_list, ~colData(.) %>% tbl_df2())
gene_meta_list = purrr::map(se_list, ~rowData(.) %>% tbl_df2())

#Construct batches list
batches = (dplyr::data_frame(a = c(1:50), b = 50) %>% 
             dplyr::mutate(batch_id = paste(a,b,sep ="_")))$batch_id
batches_list = idVectorToList(batches)

#Import leafcutter results
leafcutter_files = paste("results/acLDL/diff_splicing/differential_events_batch.leafcutter",batches,"rds", sep = ".") %>%
  as.list()
leafcutter_res = purrr::map(leafcutter_files, ~readRDS(.))

#Import ensembl_87 results
ensembl_files = paste("results/acLDL/diff_splicing/differential_events_batch.ensembl_87",batches,"rds", sep = ".") %>%
  as.list()
ensembl_res = purrr::map(ensembl_files, ~readRDS(.))

#Import leafcutter results
revised_files = paste("results/acLDL/diff_splicing/differential_events_batch.revisedAnnotation",batches,"rds", sep = ".") %>%
  as.list()
revised_res = purrr::map(revised_files, ~readRDS(.))

#Extract results
leafcutter_results = extractDTUResults(leafcutter_res, se_leafcutter)
revised_results = extractDTUResults(revised_res, se_revised)
ensembl_results = extractDTUResults(ensembl_res, se_ensembl)
results_list = list(ensembl_87 = ensembl_results, revisedAnnotation = revised_results, leafcutter = leafcutter_results)
saveRDS(results_list, "acLDL_figures/tables/DTU_genes.rds")


