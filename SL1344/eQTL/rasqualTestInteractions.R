library("devtools")
library("plyr")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                  levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Load p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_qtl_df = extractQTLsFromList(rasqual_min_pvalues, fdr_cutoff = 0.1)
joint_pairs = dplyr::select(rasqual_qtl_df, gene_id, snp_id) %>% unique() 

rasqual_selected_pvalues = readRDS("results/SL1344/eQTLs/rasqual_selected_pvalues.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Filter VCF
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, filtered_vcf$genotypes, .8)

#Naive vs IFNg
covariate_names = c("PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6", "sex_binary")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                      paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(tbl_df(filtered_pairs), combined_expression_data$cqn, combined_expression_data$sample_metadata, filtered_vcf, formula_qtl, formula_interaction)
interaction_df = postProcessInteractionPvalues(interaction_results)
saveRDS(interaction_df, "results/SL1344/eQTLs/SL1344_interaction_pvalues.rds")
interaction_df = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues.rds")
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)

#Make a Q-Q plot for the interaction p-values
qq_df = dplyr::mutate(interaction_df, p_eigen = p_nominal) %>% addExpectedPvalue()
qq_plot = ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_nominal,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")
ggsave("figures/supplementary/eQTL_interaction_Q-Q_plot.pdf", plot = qq_plot, width = 4, height = 4)

#Extract effect sizes for all gene-snp pairs from RASQUAL data
beta_list = extractAndProcessBetas(dplyr::select(interaction_hits, gene_id, snp_id), rasqual_selected_pvalues, "naive")
beta_list$beta_summaries = dplyr::mutate(beta_list$beta_summaries, max_naive_ratio = max_abs_beta/abs(naive))

#Find QTLs that appear
set.seed(42)
appear_qtls = dplyr::filter(beta_list$beta_summaries, abs(naive) <= 0.59, max_abs_beta - abs(naive) >= 0.32)
appear_betas = dplyr::semi_join(beta_list$beta_summaries[,1:6], appear_qtls, by = c("gene_id", "snp_id"))
appear_clusters = clusterBetasKmeans(appear_betas, 6) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_list$beta_df, by = c("gene_id", "snp_id"))

#Reorder clusters
cluster_reorder = data_frame(cluster_id = c(1,2,3,4,5,6), new_cluster_id = c(1,4,6,2,5,3))

#Make heatmap of effect sizes
appear_betas = appear_clusters %>% 
  dplyr::group_by(gene_id, snp_id) %>% 
  dplyr::mutate(beta_scaled = beta/max(beta)) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>%
  dplyr::left_join(cluster_reorder, by = "cluster_id") %>% #Reorder clusters
  dplyr::left_join(figureNames(), by = "condition_name") #Add figure names

#Count the number of qtls
appear_count = dplyr::select(appear_betas, gene_id, snp_id) %>% unique() %>% nrow()
ylabel = paste(appear_count, "eQTLs")
effect_size_heatmap = ggplot(appear_betas, aes(x = figure_name, y = gene_name, fill = beta_scaled)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ylab(ylabel) + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Relative effect", midpoint = 0) +
  theme_light() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.x = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))
ggsave("figures/main_figures/eQTLs_appear_kmeans_heatmap.pdf",effect_size_heatmap, width = 3.5, height = 4)
ggsave("figures/main_figures/eQTLs_appear_kmeans_heatmap.png",effect_size_heatmap, width = 3.5, height = 4)


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
set.seed(41)
disappear_qtls = dplyr::filter(beta_list$beta_summaries, abs(naive) > 0.59, abs(naive) - min_abs_beta > 0.32)
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

disappear_betas = disappear_clusters %>% dplyr::group_by(gene_id, snp_id) %>% dplyr::mutate(beta_scaled = beta/max(beta)) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>%
  dplyr::semi_join(disappear_cluster_sizes, by = "cluster_id")
disappear_count = disappear_cluster_sizes$count %>% sum()
ylabel = paste(disappear_count, "eQTLs")
dis_effect_size_heatmap = ggplot(disappear_betas, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + 
  facet_grid(cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ylab(ylabel) + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Relative effect", midpoint = 0) +
  theme_grey() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), 
        axis.title.x = element_blank(), axis.text.x=element_text(angle = 15)) +
  theme(panel.margin = unit(0.2, "lines"))
ggsave("figures/supplementary/eQTLs_disappear_kmeans_heatmap.pdf",dis_effect_size_heatmap, width = 4.5, height = 3)
ggsave("figures/supplementary/eQTLs_disappear_kmeans_heatmap.png",dis_effect_size_heatmap, width = 4.5, height = 3)


disappear_means_plot = ggplot(disappear_cluster_means, aes(x = condition_name, y = beta_mean, group = cluster_id)) + 
  geom_point() + geom_line() + facet_wrap(~cluster_id) + 
  geom_errorbar(aes(ymin = beta_mean - beta_sd, ymax = beta_mean + beta_sd, width = .2)) + 
  geom_text(aes(label = paste("n = ", count, sep = ""), x = condition_name, y = 1.5), data = disappear_cluster_sizes) +
  xlab("Condition") + 
  ylab("Log2 fold-change")
ggsave("results/SL1344/eQTLs/properties/eQTLs_disappear_cluster_means.pdf",disappear_means_plot, width = 6, height = 6)

#Export clustering results
variable_qtls = list(appear = appear_betas, disappear = disappear_betas)
saveRDS(variable_qtls, "results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")




#### ASE data ####
#Fetch ASE data from disk
exon_ranges = constructExonRanges("ENSG00000144228", "rs12621644", combined_expression_data$gene_metadata)
sample_meta = dplyr::select(combined_expression_data$sample_metadata, sample_id, condition_name, genotype_id)
ase_data = fetchGeneASEData(exon_ranges, "results/SL1344/combined_ASE_counts.sorted.txt.gz", sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)


#Make plot
ase_data = fetchGeneASEData(exon_ranges, "results/SL1344/combined_ASE_counts.sorted.txt.gz", sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)
plotting_data = filterASEforPlotting(ase_data) %>% dplyr::filter(lead_snp_value == 1)
new_data_plot = ggplot(plotting_data, aes(x = factor(lead_snp_value), y = abs(0.5-ratio), label = sample_id)) + 
  facet_grid(feature_snp_id~condition_name) +
  geom_boxplot(outlier.shape = NA) + 
  geom_text() +
  geom_jitter(position = position_jitter(width = .1)) +
  xlab("Feature SNP id") + 
  ylab("Allelic imbalance")
ggsave("results/SL1344/SPOPL_AI_new.pdf", new_data_plot, width = 8, height = 10)


new_samples = unique(b$genotype_id)
old_sample_meta = dplyr::filter(sample_meta, !(sample_meta$genotype_id %in% new_samples))

ase_data_old = fetchGeneASEData(exon_ranges, "results/SL1344/combined_ASE_counts.sorted.old.txt.gz", old_sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)
plotting_data = filterASEforPlotting(ase_data_old) %>% dplyr::filter(lead_snp_value == 1)
old_data_plot = ggplot(plotting_data, aes(x = factor(lead_snp_value), y = abs(0.5-ratio), label = sample_id)) + 
  facet_grid(feature_snp_id~condition_name) +
  geom_boxplot(outlier.shape = NA) +
  geom_text() +
  geom_jitter(position = position_jitter(width = .1)) +
  xlab("Feature SNP id") + 
  ylab("Allelic imbalance")
ggsave("results/SL1344/SPOPL_AI_old.pdf", old_data_plot, width = 8, height = 10)

SPOPL_read_counts = plotEQTL("ENSG00000144228", "rs12621644", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
ggsave("results/SL1344/SPOPL_between_individual.pdf", SPOPL_read_counts, width = 8, height = 8)


SPOPL_read_counts = plotEQTL("ENSG00000168310", "rs34156200", combined_expression_data$cqn, vcf_file$genotypes, 
                             combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)

plotEQTL("ENSG00000005844", "rs11574938", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)


#IRF2 eQTL
plotEQTL("ENSG00000168310", "rs13149699", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)



