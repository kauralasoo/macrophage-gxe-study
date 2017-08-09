library("dplyr")
library("tidyr")
library("purrr")
library("UpSetR")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import SummarizedExperiment
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")

#Identify genes in the MHC region that should be excluded
mhc_fairfax = dplyr::filter(tbl_df2(rowData(se_fairfax)), chr == "6", gene_start > 28477797, gene_end < 33448354) %>%
  dplyr::rename(phenotype_id = probe_id)

#Gene names
probe_name_map = dplyr::select(tbl_df2(rowData(se_fairfax)), probe_id, gene_name) %>% 
  dplyr::rename(phenotype_id = probe_id)

#Import GWAS traits
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt",
                                     col_names = c("trait","file_name", "type")) %>%
  dplyr::filter(!(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "CD_2012", 
                               "RA_2012", "T2D_1", "MS", "T1D", "T1D_2", "PBC", "CAD_2017"))) %>%
  dplyr::filter(type == "Autoimmune")

#Import coloc output for each subset of data
full_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".full.2e5.txt", 
                                              coloc_prefix = "processed/Fairfax/coloc/",
                                              PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                              gwas_pval_thresh = 1e-6, mhc_phenotypes = mhc_fairfax)$coloc_filtered %>%
  dplyr::left_join(probe_name_map, by = "phenotype_id") %>%
  dplyr::select(-.row)

shared_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".shared.2e5.txt", 
                                           coloc_prefix = "processed/Fairfax/coloc/",
                                           PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                           gwas_pval_thresh = 1e-6, mhc_phenotypes = mhc_fairfax)$coloc_filtered %>%
  dplyr::left_join(probe_name_map, by = "phenotype_id") %>%
  dplyr::select(-.row)

shared_84_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".shared_84.2e5.txt", 
                                           coloc_prefix = "processed/Fairfax/coloc/",
                                           PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                           gwas_pval_thresh = 1e-6, mhc_phenotypes = mhc_fairfax)$coloc_filtered %>%
  dplyr::left_join(probe_name_map, by = "phenotype_id") %>%
  dplyr::select(-.row)


#Keep one gene per summarised trait
condensed_eQTL_hits = dplyr::group_by(shared_84_200kb_hits, summarised_trait, phenotype_id) %>% 
  dplyr::arrange(summarised_trait, -PP.H4.abf) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(summarised_trait, gwas_lead, -PP.H4.abf) %>% 
  dplyr::group_by(summarised_trait, gwas_lead) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()
shared_84_200kb_filtered_hits = dplyr::semi_join(shared_84_200kb_hits, condensed_eQTL_hits, by = c("phenotype_id", "trait"))# %>%
  #dplyr::filter(gene_name != "FADS2")

#Keep one gene per summarised trait
condensed_full_eQTL_hits = dplyr::group_by(full_200kb_hits, summarised_trait, phenotype_id) %>% 
  dplyr::arrange(summarised_trait, -PP.H4.abf) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(summarised_trait, gwas_lead, -PP.H4.abf) %>% 
  dplyr::group_by(summarised_trait, gwas_lead) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()
full_200kb_filtered_hits = dplyr::semi_join(full_200kb_hits, condensed_full_eQTL_hits, by = c("phenotype_id", "trait"))# %>%
#dplyr::filter(gene_name != "FADS2")


#Map to Salmonella condition names
name_map = data_frame(condition_name = c("CD14","IFN","LPS2","LPS24"), condition_name2 = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
name_map2 = data_frame(condition_name = c("CD14","IFN","LPS2","LPS24"), condition_name2 = c("N", "I", "S", "I+S"))

#Partition into conditions
renamed_olaps = dplyr::left_join(shared_84_200kb_filtered_hits, name_map, by = "condition_name") %>%
  dplyr::mutate(condition_name = condition_name2) %>% dplyr::select(-condition_name2)
eqtl_coloc_counts = countConditionSpecificOverlaps(renamed_olaps, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
eqtl_total_counts = group_by(eqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) %>% 
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "84 samples") %>%
  dplyr::rename(condition_name2 = figure_name) %>%
  dplyr::left_join(name_map2, by = "condition_name2") %>%
  dplyr::mutate(condition_name = ifelse(condition_name == "IFN", "CD14\nIFN", condition_name)) %>%
  dplyr::mutate(condition_name = ifelse(condition_name == "LPS2", "CD14\nIFN\nLPS2", condition_name)) %>%
  dplyr::mutate(condition_name = ifelse(condition_name == "LPS24", "CD14\nIFN\nLPS2\nLPS24", condition_name))

#Count colocs in each condition separately
rna_condition_count = dplyr::filter(shared_84_200kb_filtered_hits, PP_power > 0.8, PP_coloc > 0.9) %>% 
  dplyr::select(summarised_trait, gene_name, condition_name) %>% 
  unique() %>% 
  dplyr::group_by(condition_name) %>% 
  dplyr::summarise(coloc_count = length(condition_name)) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(phenotype = "84 samples")

coloc_condition_plot = ggplot(rna_condition_count, aes(x = condition_name, y = coloc_count, group = phenotype, fill = phenotype)) + 
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Condition") + 
  ylab("Number of overlaps") +
  scale_y_continuous(limits = c(0,25)) +
  theme_light() + 
  scale_fill_manual(values = c("#e66101","#5e3c99"), name = "") +
  theme(legend.position = "top")
ggsave("figures/supplementary/fairfax_coloc_QTL_condition_counts.pdf", plot = coloc_condition_plot, width = 2.6, height = 3)


#### Count overlaps for the full dataset
#Partition into conditions
renamed_olaps = dplyr::left_join(full_200kb_filtered_hits, name_map, by = "condition_name") %>%
  dplyr::mutate(condition_name = condition_name2) %>% dplyr::select(-condition_name2)
eqtl_coloc_counts = countConditionSpecificOverlaps(renamed_olaps, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
eqtl_total_counts_full = group_by(eqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) %>% 
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "Full dataset") %>%
  dplyr::rename(condition_name2 = figure_name) %>%
  dplyr::left_join(name_map2, by = "condition_name2") %>%
  dplyr::mutate(condition_name = ifelse(condition_name == "IFN", "CD14\nIFN", condition_name)) %>%
  dplyr::mutate(condition_name = ifelse(condition_name == "LPS2", "CD14\nIFN\nLPS2", condition_name)) %>%
  dplyr::mutate(condition_name = ifelse(condition_name == "LPS24", "CD14\nIFN\nLPS2\nLPS24", condition_name))

#Count colocs in each condition separately
rna_condition_count_full = dplyr::filter(full_200kb_filtered_hits, PP_power > 0.8, PP_coloc > 0.9) %>% 
  dplyr::select(summarised_trait, gene_name, condition_name) %>% 
  unique() %>% 
  dplyr::group_by(condition_name) %>% 
  dplyr::summarise(coloc_count = length(condition_name)) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(phenotype = "Full dataset")

#Make a barplot with overlap counts
coloc_counts_plot = ggplot(bind_rows(eqtl_total_counts, eqtl_total_counts_full), aes(x = condition_name, y = total_overlap, group = phenotype, color = phenotype)) + 
  geom_point() +
  geom_line() +
  xlab("Conditions included") + 
  ylab("Cumulative number of overlaps") +
  scale_y_continuous(limits = c(0,110)) +
  theme_light() + 
  scale_color_manual(values = c("#e66101","#5e3c99"), name = "") +
  theme(legend.position = "top")
ggsave("figures/supplementary/fairfax_coloc_QTL_counts.pdf", plot = coloc_counts_plot, width = 2.6, height = 3)


coloc_condition_plot = ggplot(bind_rows(rna_condition_count, rna_condition_count_full), aes(x = condition_name, y = coloc_count, group = phenotype, fill = phenotype)) + 
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Condition") + 
  ylab("Number of overlaps") +
  scale_y_continuous(limits = c(0,25)) +
  theme_light() + 
  scale_fill_manual(values = c("#e66101","#5e3c99"), name = "") +
  theme(legend.position = "top")
ggsave("figures/supplementary/fairfax_coloc_QTL_condition_counts.pdf", plot = coloc_condition_plot, width = 2.6, height = 3)



#Estimate how many of the "stimulated" hits are detected in the naive condition with 414 samples (5x the sample size)
gained_hits = dplyr::filter(eqtl_coloc_counts, condition_name != "naive") %>%
  dplyr::left_join(probe_name_map) %>%
  dplyr::select(summarised_trait, gene_name) %>%
  unique()

#Add gwas lead variant
gained_hits_lead = dplyr::semi_join(shared_84_200kb_hits, gained_hits, by = c("summarised_trait", "gene_name")) %>%
  dplyr::select(summarised_trait, gene_name, gwas_lead) %>% unique()

full_naive_hits = dplyr::filter(full_200kb_hits, condition_name == "CD14", PP_power > 0.8, PP_coloc > 0.9)
discovered_in_full = dplyr::semi_join(full_naive_hits, gained_hits_lead, by = c("summarised_trait", "gene_name", "gwas_lead"))

#Calculate fraction
nrow(discovered_in_full)
nrow(gained_hits_lead)
nrow(discovered_in_full)/nrow(gained_hits_lead)


