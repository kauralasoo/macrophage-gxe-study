library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                  levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Load p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_min_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = ldply(rasqual_min_hits, .id = "condition_name")
joint_pairs = dplyr::select(min_pvalue_df, gene_id, snp_id) %>% unique() 

rasqual_selected_pvalues = readRDS("results/SL1344/eQTLs/rasqual_selected_pvalues.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Calculate R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, filtered_vcf$genotypes, .8)

#Naive vs IFNg
covariate_names = c("sex_binary", "ng_ul_mean","macrophage_diff_days","rna_auto", "max_purity_filtered", "harvest_stimulation_days",
                    "PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                      paste(covariate_names, collapse = " + "), sep = "+ "))
#Test for interactions
interaction_results = testMultipleInteractions(filtered_pairs, combined_expression_data, filtered_vcf, formula_qtl, formula_interaction)
interaction_df = postProcessInteractionPvalues(interaction_results)
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)

#Extract effect sizes for all gene-snp pairs from RASQUAL data
beta_list = extractAndProcessBetas(dplyr::select(interaction_hits, gene_id, snp_id), rasqual_selected_pvalues, "naive")

#Find QTLs that appear
appear_qtls = dplyr::filter(beta_list$beta_summaries, abs(naive) <= 0.59, max_abs_diff >= 0.59, max_abs_beta >= 0.59)
appear_betas = dplyr::semi_join(beta_list$beta_summaries[,1:6], appear_qtls, by = c("gene_id", "snp_id"))
appear_clusters = clusterBetasKmeans(appear_betas, 6) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_list$beta_df, by = c("gene_id", "snp_id"))
appear_plot = ggplot(appear_clusters, aes(x = condition_name, y = beta, group = paste(gene_id, snp_id))) + 
  geom_line() + facet_wrap(~cluster_id)
ggsave("results/SL1344/eQTLs/properties/eQTLs_appear_kmeans.pdf",appear_plot, width = 10, height = 10)

#Calculate mean effect size in each cluster and condition
appear_cluster_means = calculateClusterMeans(appear_clusters)
cluster_sizes = calculateClusterSizes(appear_clusters)

appear_means_plot = ggplot(appear_cluster_means, aes(x = condition_name, y = beta_mean, group = cluster_id)) + 
  geom_point() + geom_line() + facet_wrap(~cluster_id) + 
  geom_errorbar(aes(ymin = beta_mean - beta_sd, ymax = beta_mean + beta_sd, width = .2)) + 
  geom_text(aes(label = paste("n = ", count, sep = ""), x = condition_name, y = 1.5), data = cluster_sizes) +
  xlab("Condition") + 
  ylab("Log2 fold-change")
ggsave("results/SL1344/eQTLs/properties/eQTLs_appear_cluster_means.pdf",appear_means_plot, width = 9, height = 6)

#Find QTLs that disappear
#Look for QTLs that disappear after stimulation
disappear_qtls = dplyr::filter(beta_list$beta_summaries, abs(naive) > 0.59, max_abs_diff >= 0.59, min_abs_beta <= 0.59)
disappear_betas = dplyr::semi_join(beta_list$beta_summaries[,1:6], disappear_qtls, by = c("gene_id", "snp_id"))
disappear_clusters = clusterBetasKmeans(disappear_betas, 7) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_list$beta_df, by = c("gene_id", "snp_id"))
disappear_plot = ggplot(disappear_clusters, aes(x = condition_name, y = abs(beta), group = paste(gene_id, snp_id))) + 
  geom_line() + facet_wrap(~cluster_id)
ggsave("results/SL1344/eQTLs/properties/eQTLs_disappear_kmeans.pdf",disappear_plot, width = 10, height = 10)

#Calculate mean effect size in each cluster and condition
disappear_cluster_sizes = calculateClusterSizes(disappear_clusters, selected_condition = "IFNg_SL1344") %>% 
  dplyr::filter(count > 10)
disappear_cluster_means = calculateClusterMeans(disappear_clusters) %>%
  dplyr::semi_join(disappear_cluster_sizes, by = "cluster_id")

disappear_means_plot = ggplot(disappear_cluster_means, aes(x = condition_name, y = beta_mean, group = cluster_id)) + 
  geom_point() + geom_line() + facet_wrap(~cluster_id) + 
  geom_errorbar(aes(ymin = beta_mean - beta_sd, ymax = beta_mean + beta_sd, width = .2)) + 
  geom_text(aes(label = paste("n = ", count, sep = ""), x = condition_name, y = 1.5), data = disappear_cluster_sizes) +
  xlab("Condition") + 
  ylab("Log2 fold-change")
ggsave("results/SL1344/eQTLs/properties/eQTLs_disappear_cluster_means.pdf",disappear_means_plot, width = 6, height = 6)




#### Find most associated peaks for each gene ####
rna_betas = extractAndProcessBetas(dplyr::select(interaction_hits, gene_id, snp_id), rasqual_selected_pvalues, "naive")
rna_appear_qtls = dplyr::filter(rna_betas$beta_summaries, abs(naive) <= 0.59, max_abs_diff >= 0.59, max_abs_beta >= 0.59)

#Extract most associated peaks for each of the eQTL genes
appear_hits = dplyr::left_join(appear_qtls, gene_name_map, by = "gene_id")

#Construct GRanges object of SNP positions
selected_snps = dplyr::semi_join(filtered_vcf$snpspos, rna_appear_qtls, by = c("snpid" = "snp_id")) %>%
  dplyr::transmute(snp_id = snpid, seqnames = chr, start = pos, end = pos, strand = "*") %>%
  dataFrameToGRanges()

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Fetch corresponding SNPs from ATAC data
atac_tabix_list = list(naive = "~/projects/macrophage-chromatin/naive_100kb.sorted.unskipped.txt.gz",
                       IFNg = "~/projects/macrophage-chromatin/IFNg_100kb.sorted.unskipped.txt.gz",
                       SL1344 = "~/projects/macrophage-chromatin/SL1344_100kb.sorted.unskipped.txt.gz",
                       IFNg_SL1344 = "~/projects/macrophage-chromatin/IFNg_SL1344_100kb.sorted.unskipped.txt.gz")
atac_rasqual_list = lapply(atac_tabix_list, function(tabix, snps) tabixFetchSNPs(snps, tabix), selected_snps)

#Identify QTLs that appeat after specific stimuli
ifng_appear_qtls = dplyr::filter(rna_appear_qtls, abs(IFNg) >= 0.59, abs(IFNg_diff) >= 0.59) %>%
  dplyr::group_by(gene_id) %>% dplyr::arrange(-max_abs_beta) %>%
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()
sl1344_appear_qtls = dplyr::filter(rna_appear_qtls, abs(SL1344) >= 0.59, abs(SL1344_diff) >= 0.59) %>%
  dplyr::group_by(gene_id) %>% dplyr::arrange(-max_abs_beta) %>%
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()
ifng_sl1344_appear_qtls = dplyr::filter(appear_clusters, cluster_id == 6) %>% 
  dplyr::select(gene_id, snp_id) %>% unique() %>% 
  dplyr::left_join(rna_appear_qtls, by = c("gene_id","snp_id"))

#Find candidate ATAC peaks that overlap with condition-specific eQTLs
ifng_peaks = findMostAssociatedPeakPerSNP(ifng_appear_qtls, atac_snp_tables[["IFNg"]]) %>% dplyr::filter(p_fdr < 0.1)
sl1344_peaks = findMostAssociatedPeakPerSNP(sl1344_appear_qtls, atac_snp_tables$SL1344) %>% dplyr::filter(p_fdr < 0.1)
ifng_sl1344_peaks = findMostAssociatedPeakPerSNP(ifng_sl1344_appear_qtls, atac_snp_tables$IFNg_SL1344) %>% dplyr::filter(p_fdr < 0.1)

#IFNg
ifng_effects = prepareBetasDf(ifng_appear_qtls, rna_betas, atac_rasqual_list, gene_name_map, 
                              appear_condition = "IFNg", rank_by = "IFNg_diff") %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  dplyr::group_by(snp_id, peak_id, gene_id, phenotype) %>% 
  dplyr::mutate(beta_std = (beta - mean(beta))/sd(beta), beta_scaled = beta/max(beta)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_binary = ifelse(beta >= 0.59, 1, 0))

effect_size_heatmap = ggplot(ifng_effects, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + facet_wrap(~phenotype) + geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = scales::muted("red"), mid = "white", high = scales::muted("blue"), name = "Beta", midpoint = 0) 
ggsave("results/SL1344/eQTLs/properties/eQTLs_vs_caQTL_IFNg_heatmap.pdf", effect_size_heatmap, width = 5, height = 7)

#SL1344
sl1344_effects = prepareBetasDf(sl1344_appear_qtls, rna_betas, atac_rasqual_list, gene_name_map, 
                              appear_condition = "SL1344", rank_by = "SL1344_diff") %>%
  dplyr::filter(condition_name %in% c("naive","SL1344")) %>%
  dplyr::group_by(snp_id, peak_id, gene_id, phenotype) %>% 
  dplyr::mutate(beta_std = (beta - mean(beta))/sd(beta), beta_scaled = beta/max(beta)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_binary = ifelse(beta >= 0.59, 1, 0))

effect_size_heatmap = ggplot(sl1344_effects, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + facet_wrap(~phenotype) + geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = scales::muted("red"), mid = "white", high = scales::muted("blue"), name = "Beta", midpoint = 0) 
ggsave("results/SL1344/eQTLs/properties/eQTLs_vs_caQTL_SL1344_heatmap.pdf", effect_size_heatmap, width = 5, height = 7)

#SL1344
ifng_sl1344_effects = prepareBetasDf(ifng_sl1344_appear_qtls, rna_betas, atac_rasqual_list, gene_name_map, 
                                appear_condition = "IFNg_SL1344", rank_by = "IFNg_SL1344_diff") %>%
  dplyr::group_by(snp_id, peak_id, gene_id, phenotype) %>% 
  dplyr::mutate(beta_std = (beta - mean(beta))/sd(beta), beta_scaled = beta/max(beta)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_binary = ifelse(beta >= 0.59, 1, 0))

effect_size_heatmap = ggplot(ifng_sl1344_effects, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + facet_wrap(~phenotype) + geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = scales::muted("red"), mid = "white", high = scales::muted("blue"), name = "Beta", midpoint = 0) 
ggsave("results/SL1344/eQTLs/properties/eQTLs_vs_caQTL_IFNg_SL1344_heatmap.pdf", effect_size_heatmap, width = 5, height = 4)


#Make plots
makeMultiplePlots(ifng_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_IFNg/", 7,7)
makeMultiplePlots(sl1344_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_SL1344/", 7,7)
makeMultiplePlots(ifng_sl1344_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_IFNg_SL1344/", 7,7)


#Fetch ASE data from disk
exon_ranges = constructExonRanges("ENSG00000144228", "rs12621644", combined_expression_data$gene_metadata)
sample_meta = dplyr::select(combined_expression_data$sample_metadata, sample_id, condition_name, genotype_id)
ase_data = fetchGeneASEData(exon_ranges, "results/SL1344/combined_ASE_counts.sorted.txt.gz", sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)


#Make plot
plotting_data = filterASEforPlotting(ase_data) %>% dplyr::filter(total_count > 10) %>% dplyr::filter(lead_snp_value == 1)
ggplot(plotting_data, aes(x = factor(lead_snp_value), y = abs(0.5-ratio))) + 
  facet_grid(feature_snp_id~condition_name) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = .1)) +
  xlab("Feature SNP id") + 
  ylab("Reference allele ratio")


