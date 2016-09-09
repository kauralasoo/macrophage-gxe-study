library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import event dataset
se_ensembl = readRDS("results/SL1344/combined_reviseAnnotations_transcript_quants.rds")
event_metadata = rowData(se_ensembl)

#Extract event names
event_names = dplyr::select(tbl_df(as.data.frame(event_metadata)), 
                            gene_id, ensembl_gene_id, event_type) %>% unique()

#Import differentially expressed events
event_columns = c("gene_id","transcript_id","naive","IFNg","SL1344","IFNg_SL1344",
                  "SL1344_lr","SL1344_pvalue","IFNg_lr","IFNg_pvalue","IFNg_SL1344_lr",
                  "IFNg_SL1344_pvalue","all_lr","all_pvalue")
diff_events = readr::read_delim("results/SL1344/diff_splicing/differential_events.txt.gz", 
                                col_names = event_columns, delim ="\t", col_types = "ccdddddddddddd")

#Calculate differences in proportion
diff_events = dplyr::mutate(diff_events, SL1344_diff = naive - SL1344, 
                            IFNg_diff = naive - IFNg, IFNg_SL1344_diff = naive - IFNg_SL1344)

#Salmonella DE transcripts
SL1344_events = dplyr::filter(diff_events, SL1344_pvalue < 0.01) %>% 
  dplyr::left_join(event_names, by = "gene_id") %>%
  dplyr::select(transcript_id,gene_id, ensembl_gene_id, event_type, SL1344_diff, SL1344_pvalue) %>% 
  dplyr::mutate(SL1344_abs_diff = abs(SL1344_diff)) %>% 
  dplyr::group_by(ensembl_gene_id, event_type) %>% 
  dplyr::arrange(-SL1344_abs_diff) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() 

sl1344_hits = dplyr::filter(SL1344_events, SL1344_abs_diff > 0.1) #At least 10% change in proportions

#IFNg DE events
IFNg_events = dplyr::filter(diff_events, IFNg_pvalue < 0.01) %>% 
  dplyr::left_join(event_names, by = "gene_id") %>%
  dplyr::select(transcript_id,gene_id, ensembl_gene_id, event_type, IFNg_diff, IFNg_pvalue) %>% 
  dplyr::mutate(IFNg_abs_diff = abs(IFNg_diff)) %>% 
  dplyr::group_by(ensembl_gene_id, event_type) %>% 
  dplyr::arrange(-IFNg_abs_diff) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() 

IFNg_hits = dplyr::filter(IFNg_events, IFNg_abs_diff > 0.1) #At least 10% change in proportions

#IFNg_SL1344 DE events
IFNg_SL1344_events = dplyr::filter(diff_events, IFNg_SL1344_pvalue < 0.01) %>% 
  dplyr::left_join(event_names, by = "gene_id") %>%
  dplyr::select(transcript_id,gene_id, ensembl_gene_id, event_type, IFNg_SL1344_diff, IFNg_SL1344_pvalue) %>% 
  dplyr::mutate(IFNg_SL1344_abs_diff = abs(IFNg_SL1344_diff)) %>% 
  dplyr::group_by(ensembl_gene_id, event_type) %>% 
  dplyr::arrange(-IFNg_SL1344_abs_diff) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() 

IFNg_SL1344_hits = dplyr::filter(IFNg_SL1344_events, IFNg_SL1344_abs_diff > 0.1) #At least 10% change in proportions

#Make some plots of effect sizes
ggplot(SL1344_events, aes(x = event_type, y = SL1344_abs_diff)) + geom_violin() + geom_boxplot(width = 0.1)
ggplot(IFNg_events, aes(x = event_type, y = IFNg_abs_diff)) + geom_violin() + geom_boxplot(width = 0.1)
ggplot(IFNg_SL1344_events, aes(x = event_type, y = IFNg_SL1344_abs_diff)) + geom_violin() + geom_boxplot(width = 0.1)

