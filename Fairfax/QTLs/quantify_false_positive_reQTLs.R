library("readr")
library("dplyr")
library("SummarizedExperiment")

importFairfaxAFC <- function(dir, suffix = ".aFC_results.txt"){
  #Import aFC estimates in each condition for the full dataset (228 ind)
  cd14_table = readr::read_tsv(file.path(dir, paste0("CD14", suffix))) %>% 
    dplyr::transmute(gene_id = pid, snp_id = sid, CD14 = log2_aFC)
  ifn_table = readr::read_tsv(file.path(dir, paste0("IFN", suffix))) %>% 
    dplyr::transmute(gene_id = pid, snp_id = sid, IFN = log2_aFC)
  lps2_table = readr::read_tsv(file.path(dir, paste0("LPS2", suffix))) %>% 
    dplyr::transmute(gene_id = pid, snp_id = sid, LPS2 = log2_aFC)
  lps24_table = readr::read_tsv(file.path(dir, paste0("LPS24", suffix))) %>% 
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
gene_name_map = rowData(se_fairfax) %>% tbl_df2() %>% 
  dplyr::mutate(gene_id = probe_id) %>%
  dplyr::select(gene_id, gene_name)

#Import results from interaction test
shared_interactions = readRDS("results/Fairfax/interactions/shared_interactions_lme4.rds")
shared_hits = dplyr::filter(shared_interactions, p_fdr < 0.1)

#Import aFCs
shared_diff_df = importFairfaxAFC("processed/Fairfax/qtltools/input/shared/")
saveRDS(shared_diff_df, "results/Fairfax/interactions/shared_aFC_estimates.rds")

#Import aFC values for subset QTLs
shared_diff_df_84 = importFairfaxAFC("processed/Fairfax/qtltools/input/shared/", suffix = ".shared_84.aFC_results.txt")
shared_diff_df_42 = importFairfaxAFC("processed/Fairfax/qtltools/input/shared/", suffix = ".shared_42.aFC_results.txt")

#Extract the correct columns
shared_effects_84 = dplyr::select(shared_diff_df_84, gene_id, snp_id, CD14, IFN, LPS2, LPS24) %>% 
  tidyr::gather("condition_name", "aFC", CD14:LPS24)
shared_effects_42 = dplyr::select(shared_diff_df_42, gene_id, snp_id, CD14, IFN, LPS2, LPS24) %>% 
  tidyr::gather("condition_name", "aFC", CD14:LPS24)

#Identify condition-specific QTLs
shared_fc_hits = dplyr::semi_join(shared_diff_df, shared_hits, by = c("gene_id", "snp_id")) %>% 
  dplyr::filter(abs(CD14) <= 0.29, max_fc >= 0.29, max_diff >= 0.29)
shared_fc_hits_20 = dplyr::semi_join(shared_diff_df, shared_hits, by = c("gene_id", "snp_id")) %>% 
  dplyr::filter(abs(CD14) <= 0.29*1.2, max_fc >= 0.29*0.8, max_diff >= 0.29*0.8)


##### Import subset interaction results #####
shared_84_interactions = readRDS("results/Fairfax/interactions/shared_84_interactions_lme4.rds")
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

#Compare fold_change distributions
subset_effects = dplyr::select(subset_fc_hits, gene_id, snp_id, CD14, IFN, LPS2, LPS24) %>% 
  tidyr::gather("condition_name", "aFC", CD14:LPS24)
max_effects = dplyr::group_by(subset_effects, gene_id, snp_id) %>% 
  dplyr::mutate(aFC_abs = abs(aFC)) %>% 
  dplyr::arrange(gene_id, snp_id, -aFC_abs) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::mutate(max_condition = "stimulated (max effect)") %>% dplyr::ungroup()
naive_effects = dplyr::filter(subset_effects, condition_name == "CD14") %>%
  dplyr::mutate(aFC_abs = abs(aFC), max_condition = "naive")
effects = dplyr::bind_rows(max_effects, naive_effects)
original_effects = ggplot(effects, aes(x = aFC_abs*2)) + 
  geom_histogram() + 
  facet_wrap(~max_condition) +
  geom_vline(xintercept = 0.58, color = "red") +
  theme_light() +
  xlab("Absolute fold change") +
  ylab("Number of eQTLs")
ggsave("figures/supplementary/fairfax_eQTL_84_sample_effects.pdf", plot = original_effects, width = 4.5, height = 2.5)

#Extract the same probes from the large dataset:
shared_selected_effects = dplyr::semi_join(shared_effects_84, effects, by = c("gene_id", "snp_id", "condition_name")) %>% 
  dplyr::arrange(gene_id) %>% 
  dplyr::mutate(aFC_abs = abs(aFC), max_condition = ifelse(condition_name == "CD14", "naive", "stimulated (max effect)"))
large_sample_effects = ggplot(shared_selected_effects, aes(x = aFC_abs*2)) + 
  geom_histogram() + 
  facet_wrap(~max_condition) +
  geom_vline(xintercept = 0.58, color = "red") +
  xlab("Absolute fold change") +
  ylab("Number of eQTLs") + 
  theme_light()
ggsave("figures/supplementary/fairfax_eQTL_84_228_sample_effects.pdf", plot = large_sample_effects, width = 4.5, height = 2.5)


##### Import subset interaction results (n = 42) #####
shared_42_interactions = readRDS("results/Fairfax/interactions/shared_42_interactions_lme4.rds")
subset_hits = dplyr::filter(shared_42_interactions, p_fdr < 0.1)

#Import aFCs
subset_diff_df = importFairfaxAFC("processed/Fairfax/qtltools/input/shared_42/") 

#Filter by fold change
subset_fc_hits = dplyr::semi_join(subset_diff_df, subset_hits, by = c("gene_id", "snp_id")) %>% 
  dplyr::filter(abs(CD14) <= 0.29, max_fc >= 0.29, max_diff >= 0.29)

#Estimate the rate of false positives
length(intersect(subset_hits$gene_id, shared_hits$gene_id))/length(unique(subset_hits$gene_id))
length(intersect(subset_fc_hits$gene_id, shared_fc_hits$gene_id))/length(unique(subset_fc_hits$gene_id))
length(intersect(subset_fc_hits$gene_id, shared_fc_hits_20$gene_id))/length(unique(subset_fc_hits$gene_id))

#Compare fold_change distributions
subset_effects = dplyr::select(subset_fc_hits, gene_id, snp_id, CD14, IFN, LPS2, LPS24) %>% 
  tidyr::gather("condition_name", "aFC", CD14:LPS24)
max_effects = dplyr::group_by(subset_effects, gene_id, snp_id) %>% 
  dplyr::mutate(aFC_abs = abs(aFC)) %>% 
  dplyr::arrange(gene_id, snp_id, -aFC_abs) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::mutate(max_condition = "stimulated (max effect)") %>% dplyr::ungroup()
naive_effects = dplyr::filter(subset_effects, condition_name == "CD14") %>%
  dplyr::mutate(aFC_abs = abs(aFC), max_condition = "naive")
effects = dplyr::bind_rows(max_effects, naive_effects)
original_effects = ggplot(effects, aes(x = aFC_abs*2)) + 
  geom_histogram() + 
  facet_wrap(~max_condition) +
  geom_vline(xintercept = 0.58, color = "red") +
  theme_light() +
  xlab("Absolute fold change") +
  ylab("Number of eQTLs")
ggsave("figures/supplementary/fairfax_eQTL_42_sample_effects.pdf", plot = original_effects, width = 4.5, height = 2.5)

#Extract the same probes from the large dataset:
shared_selected_effects = dplyr::semi_join(shared_effects_42, effects, by = c("gene_id", "snp_id", "condition_name")) %>% 
  dplyr::arrange(gene_id) %>% 
  dplyr::mutate(aFC_abs = abs(aFC), max_condition = ifelse(condition_name == "CD14", "naive", "stimulated (max effect)"))
large_sample_effects = ggplot(shared_selected_effects, aes(x = aFC_abs*2)) + 
  geom_histogram() + 
  facet_wrap(~max_condition) +
  geom_vline(xintercept = 0.58, color = "red") +
  xlab("Absolute fold change") +
  ylab("Number of eQTLs") + 
  theme_light()
ggsave("figures/supplementary/fairfax_eQTL_42_228_sample_effects.pdf", plot = large_sample_effects, width = 4.5, height = 2.5)

#Export filtered interaction results
shared_final_hits = dplyr::left_join(shared_fc_hits, gene_name_map, by = "gene_id")
saveRDS(shared_final_hits, "results/Fairfax/interactions/shared_interaction_hits.rds")


