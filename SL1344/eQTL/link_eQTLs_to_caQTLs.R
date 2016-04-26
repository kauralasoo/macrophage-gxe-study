library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import RASQUAL results for the eQTLs
rasqual_selected_pvalues = readRDS("results/SL1344/eQTLs/rasqual_selected_pvalues.rds")

#Import interaction eQTLs
interaction_df = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues.rds")
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)
qtl_clusters = readRDS("results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#### Find most associated peaks for each gene ####
rna_betas = extractAndProcessBetas(dplyr::select(interaction_hits, gene_id, snp_id), rasqual_selected_pvalues, "naive")
rna_appear_qtls = dplyr::filter(rna_betas$beta_summaries, abs(naive) <= 0.59, max_abs_beta - abs(naive) >= 0.32)

#Extract most associated peaks for each of the eQTL genes
appear_hits = dplyr::left_join(rna_appear_qtls, gene_name_map, by = "gene_id")

#Construct GRanges object of SNP positions
selected_snps = dplyr::filter(vcf_file$snpspos, snpid %in% rna_appear_qtls$snp_id) %>%
  dplyr::transmute(snp_id = snpid, seqnames = chr, start = pos, end = pos, strand = "*") %>%
  dataFrameToGRanges()

#Fetch corresponding SNPs from ATAC data
atac_tabix_list = list(naive = "../macrophage-chromatin/results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz",
                       IFNg = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz",
                       SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz",
                       IFNg_SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz")
atac_snp_tables = lapply(atac_tabix_list, function(tabix, snps) rasqualTools::tabixFetchSNPs(snps, tabix), selected_snps)

#Identify QTLs that appeat after specific stimuli
ifng_appear_qtls = dplyr::filter(rna_appear_qtls, abs(IFNg) >= 0.32, abs(IFNg_diff) >= 0.32) %>%
  dplyr::group_by(gene_id) %>% dplyr::arrange(-max_abs_beta) %>%
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()
sl1344_appear_qtls = dplyr::filter(rna_appear_qtls, abs(SL1344) >= 0.32, abs(SL1344_diff) >= 0.32) %>%
  dplyr::group_by(gene_id) %>% dplyr::arrange(-max_abs_beta) %>%
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()
ifng_sl1344_appear_qtls = dplyr::filter(qtl_clusters$appear, cluster_id == 3) %>% 
  dplyr::select(gene_id, snp_id) %>% ungroup() %>% unique() %>% 
  dplyr::left_join(rna_appear_qtls, by = c("gene_id","snp_id"))

#IFNg - find corresponding ATAC peaks
ifng_effects = prepareBetasDf(ifng_appear_qtls, rna_betas, atac_snp_tables, gene_name_map, 
                              appear_condition = "IFNg", rank_by = "IFNg_diff") %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  dplyr::group_by(snp_id, peak_id, gene_id, phenotype) %>% 
  dplyr::mutate(beta_std = (beta - mean(beta))/sd(beta), beta_scaled = beta/max(beta)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_binary = ifelse(beta >= 0.59, 1, 0))

#Calculate scaled ATAC diff
scaled_diff = dplyr::filter(ifng_effects, phenotype == "ATAC") %>% 
  group_by(gene_id) %>% 
  dplyr::arrange(condition_name) %>% 
  dplyr::mutate(scaled_diff = beta_scaled[2] - beta_scaled[1]) %>% 
  dplyr::filter(condition_name == "naive") %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(scaled_diff)
ifng_effects_sorted = dplyr::mutate(ifng_effects, gene_name = factor(as.character(gene_name), levels = as.character(scaled_diff$gene_name)))

effect_size_heatmap = ggplot(ifng_effects_sorted, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + 
  facet_wrap(~phenotype) + geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Beta", midpoint = 0) 
ggsave("results/SL1344/eQTLs/properties/eQTLs_vs_caQTL_IFNg_heatmap.pdf", effect_size_heatmap, width = 5, height = 11)

#SL1344 - find corresponding ATAC peaks
sl1344_effects = prepareBetasDf(sl1344_appear_qtls, rna_betas, atac_snp_tables, gene_name_map, 
                                appear_condition = "SL1344", rank_by = "SL1344_diff") %>%
  dplyr::filter(condition_name %in% c("naive","SL1344")) %>%
  dplyr::group_by(snp_id, peak_id, gene_id, phenotype) %>% 
  dplyr::mutate(beta_std = (beta - mean(beta))/sd(beta), beta_scaled = beta/max(beta)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_binary = ifelse(beta >= 0.59, 1, 0))

#Calculate scaled ATAC diff
scaled_diff = dplyr::filter(sl1344_effects, phenotype == "ATAC") %>% 
  group_by(gene_id) %>% 
  dplyr::arrange(condition_name) %>% 
  dplyr::mutate(scaled_diff = beta_scaled[2] - beta_scaled[1]) %>% 
  dplyr::filter(condition_name == "naive") %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(scaled_diff)
sl1344_effects_sorted = dplyr::mutate(sl1344_effects, gene_name = factor(as.character(gene_name), levels = as.character(scaled_diff$gene_name)))

effect_size_heatmap = ggplot(sl1344_effects_sorted, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + facet_wrap(~phenotype) + geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Beta", midpoint = 0) 
ggsave("results/SL1344/eQTLs/properties/eQTLs_vs_caQTL_SL1344_heatmap.pdf", effect_size_heatmap, width = 5, height = 11)

#IFNg_SL1344 - find corresponding ATAC peaks
ifng_sl1344_effects = prepareBetasDf(ifng_sl1344_appear_qtls, rna_betas, atac_snp_tables, gene_name_map, 
                                     appear_condition = "IFNg_SL1344", rank_by = "IFNg_SL1344_diff") %>%
  dplyr::group_by(snp_id, peak_id, gene_id, phenotype) %>% 
  dplyr::mutate(beta_std = (beta - mean(beta))/sd(beta), beta_scaled = beta/max(beta)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_binary = ifelse(beta >= 0.59, 1, 0))

effect_size_heatmap = ggplot(ifng_sl1344_effects, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + facet_wrap(~phenotype) + geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Beta", midpoint = 0) 
ggsave("results/SL1344/eQTLs/properties/eQTLs_vs_caQTL_IFNg_SL1344_heatmap.pdf", effect_size_heatmap, width = 8, height = 7)

#Count the number of peaks per SNP
ifng_peaks_ifng = findAllAssociatedPeaksPerSNP(ifng_appear_qtls, atac_snp_tables[["IFNg"]]) %>% 
  dplyr::filter(p_bonferroni < 0.1) %>% dplyr::group_by(snp_id) %>% dplyr::summarise(ifng_peak_count = length(gene_id))
ifng_peaks_naive = findAllAssociatedPeaksPerSNP(ifng_appear_qtls, atac_snp_tables[["naive"]]) %>% 
  dplyr::filter(p_bonferroni < 0.1) %>% dplyr::group_by(snp_id) %>% dplyr::summarise(naive_peak_count = length(gene_id))
peak_count_comparison = dplyr::left_join(ifng_peaks_ifng, ifng_peaks_naive, by = "snp_id") %>% 
  dplyr::mutate(naive_peak_count = ifelse(is.na(naive_peak_count),0, naive_peak_count))
peak_count_diff = dplyr::mutate(peak_count_comparison, diff = ifng_peak_count - naive_peak_count) %>% dplyr::arrange(-diff)
write.table(peak_count_diff, "results/SL1344/eQTLs/properties/naive_ifng_ATAC_peak_count_diff.txt")


ifng_peaks_ifng = findAllAssociatedPeaksPerSNP(sl1344_appear_qtls, atac_snp_tables[["SL1344"]]) %>% 
  dplyr::filter(p_bonferroni < 0.1) %>% dplyr::group_by(snp_id) %>% dplyr::summarise(sl1344_peak_count = length(gene_id))
ifng_peaks_naive = findAllAssociatedPeaksPerSNP(sl1344_appear_qtls, atac_snp_tables[["naive"]]) %>% 
  dplyr::filter(p_bonferroni < 0.1) %>% dplyr::group_by(snp_id) %>% dplyr::summarise(naive_peak_count = length(gene_id))
peak_count_comparison = dplyr::left_join(ifng_peaks_ifng, ifng_peaks_naive, by = "snp_id") %>% 
  dplyr::mutate(naive_peak_count = ifelse(is.na(naive_peak_count),0, naive_peak_count))
peak_count_diff = dplyr::mutate(peak_count_comparison, diff = sl1344_peak_count - naive_peak_count) %>% dplyr::arrange(-diff)
write.table(peak_count_diff, "results/SL1344/eQTLs/properties/naive_ifng_ATAC_peak_count_diff.txt")

#Export gene-peak pairs
pairs = list(IFNg = ifng_effects, SL1344 = sl1344_effects, IFNg_SL1344 = ifng_sl1344_effects)
saveRDS(pairs, "results/SL1344/eQTLs/interaction_peak_gene_pairs.rds")

#CIITA eQTL
plotEQTL("ENSG00000144228", "rs145753116", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
