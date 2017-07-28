library("readr")
library("dplyr")

importFairfaxAFC <- function(dir){
  #Import aFC estimates in each condition for the full dataset (228 ind)
  cd14_table = readr::read_tsv(file.path(dir, "CD14.aFC_results.txt")) %>% 
    dplyr::transmute(gene_id = pid, snp_id = sid, CD14 = log2_aFC)
  ifn_table = readr::read_tsv(file.path(dir, "IFN.aFC_results.txt")) %>% 
    dplyr::transmute(gene_id = pid, snp_id = sid, IFN = log2_aFC)
  lps2_table = readr::read_tsv(file.path(dir, "LPS2.aFC_results.txt")) %>% 
    dplyr::transmute(gene_id = pid, snp_id = sid, LPS2 = log2_aFC)
  lps24_table = readr::read_tsv(file.path(dir, "LPS24.aFC_results.txt")) %>% 
    dplyr::transmute(gene_id = pid, snp_id = sid, LPS24 = log2_aFC)
  aFC_df = dplyr::left_join(cd14_table, ifn_table, by = c("gene_id", "snp_id")) %>%
    dplyr::left_join(lps2_table, by = c("gene_id", "snp_id")) %>%
    dplyr::left_join(lps24_table, by = c("gene_id", "snp_id"))
  
  #Calculate diffs
  shared_diff_df = dplyr::mutate(aFC_df, IFN_diff = IFN-CD14, LPS2_diff = LPS2-CD14, LPS24_diff = LPS24-CD14) %>%
    dplyr::mutate(max_fc = pmax(abs(CD14), abs(IFN), abs(LPS2), abs(LPS24))) %>%
    dplyr::mutate(max_diff = pmax(abs(LPS2_diff), abs(LPS24_diff), abs(IFN_diff)))
  
  return(shared_diff_df)
}

#Import SummarizedExperiment
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")
gene_name_map = rowData(se_shared) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id) %>%
  dplyr::select(gene_id, gene_name)

#Import results from interaction test
shared_interactions = readRDS("results/Fairfax/shared_interactions_lme4.rds")
shared_hits = dplyr::filter(shared_interactions, p_fdr < 0.1)

#Import aFCs
shared_diff_df = importFairfaxAFC("processed/Fairfax/qtltools/input/shared/")
saveRDS(shared_diff_df, "results/Fairfax/interactions/shared_aFC_estimates.rds")


#Identify condition-specific QTLs
shared_fc_hits = dplyr::semi_join(shared_diff_df, shared_hits, by = c("gene_id", "snp_id")) %>% 
  dplyr::filter(abs(CD14) <= 0.29, max_fc >= 0.29, max_diff >= 0.29)
shared_fc_hits_20 = dplyr::semi_join(shared_diff_df, shared_hits, by = c("gene_id", "snp_id")) %>% 
  dplyr::filter(abs(CD14) <= 0.29*1.2, max_fc >= 0.29*0.8, max_diff >= 0.29*0.8)


##### Import subset interaction results #####
shared_84_interactions = readRDS("results/Fairfax/shared_84_interactions_lme4.rds")
subset_hits = dplyr::filter(shared_84_interactions, p_fdr < 0.1)

#Import aFCs
subset_diff_df = importFairfaxAFC("processed/Fairfax/qtltools/input/shared_84/")

#Filter by fold change
subset_fc_hits = dplyr::semi_join(subset_diff_df, subset_hits, by = c("gene_id", "snp_id")) %>% 
  dplyr::filter(abs(CD14) <= 0.29, max_fc >= 0.29, max_diff >= 0.29)

#Estimate the rate of false positives
length(intersect(subset_hits$gene_id, shared_hits$gene_id))/length(unique(subset_hits$gene_id))
length(intersect(subset_fc_hits$gene_id, shared_fc_hits$gene_id))/length(unique(subset_fc_hits$gene_id))
length(intersect(subset_fc_hits$gene_id, shared_fc_hits_20$gene_id))/length(unique(subset_fc_hits$gene_id))

#Export filtered interaction results
shared_final_hits = dplyr::left_join(shared_fc_hits, gene_name_map, by = "gene_id")
saveRDS(shared_final_hits, "results/Fairfax/interactions/shared_interaction_hits.rds")


