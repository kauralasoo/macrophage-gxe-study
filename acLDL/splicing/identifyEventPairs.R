library("dplyr")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")

calculateMeanTranscriptExpression <- function(se){
  exp_mat = assays(se)$tpm_ratios
  mean_vector = rowMeans(exp_mat, na.rm = TRUE)
  mean_df = dplyr::data_frame(transcript_id = names(mean_vector), mean_exp = mean_vector)
  return(mean_df)
}

findOtherTranscriptId <- function(qtl_df, mean_expression_df){
  
  #Identify QTL lead transcripts
  qtl_lead_df = dplyr::transmute(qtl_df, gene_id = group_id, transcript_id = phenotype_id, qtl_lead = 1)
  
  #Identify the second most highly expressed transcript
  paired_event = dplyr::left_join(mean_expression_df, qtl_lead_df, by = c("gene_id", "transcript_id")) %>% 
    dplyr::mutate(qtl_lead = ifelse(is.na(qtl_lead), 0,qtl_lead)) %>% 
    dplyr::group_by(gene_id) %>% dplyr::mutate(has_qtl = max(qtl_lead)) %>% 
    dplyr::filter(has_qtl == 1) %>% 
    dplyr::arrange(gene_id, -qtl_lead, -mean_exp) %>% 
    dplyr::filter(row_number() <= 2) %>%
    dplyr::ungroup()
  other_phenotype_df = dplyr::filter(paired_event, qtl_lead == 0) %>% 
    dplyr::transmute(group_id = gene_id, other_phenotype_id = transcript_id)
  
  #Add the second transcript id to the QTL results
  qtl_df_other = dplyr::left_join(qtl_df, other_phenotype_df, by = "group_id")
}

meanTxExpWrapper <- function(se){
  condition_list = idVectorToList(c("Ctrl","AcLDL"))
  event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,se))
  event_metadata = tbl_df2(rowData(se)) %>%
    dplyr::select(transcript_id, gene_id)
  mean_expression = purrr::map(event_conditions, ~calculateMeanTranscriptExpression(.) %>%
                                 dplyr::left_join(event_metadata, by = "transcript_id"))
  mean_expression$Diff = mean_expression$AcLDL
  
  return(mean_expression)
}

#Import trQTL results
trQTL_min_pvalue_list = readRDS("results/acLDL/trQTLs/trQTL_min_pvalues.rds")

#Import SummarizedExperiments for all phenotypes
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
se_revised = readRDS("results/acLDL/acLDL_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/acLDL/acLDL_leafcutter_counts.rds")

#Calculate mean expression in each condition
ensembl_mean_exp = meanTxExpWrapper(se_ensembl)
revised_mean_exp = meanTxExpWrapper(se_revised)
leafcutter_mean_exp = meanTxExpWrapper(se_leafcutter)
mean_exp_list = list(ensembl_87 = ensembl_mean_exp, revisedAnnotation = revised_mean_exp, leafcutter = leafcutter_mean_exp)

#Find second most higly expressed transcript for each QTL
other_pvalues_list = purrr::map2(trQTL_min_pvalue_list, mean_exp_list, ~purrr::map2(.x, .y, ~findOtherTranscriptId(.x,.y)))
saveRDS(other_pvalues_list, "results/acLDL/trQTLs/trQTL_min_pvalues.other_tx.rds")



