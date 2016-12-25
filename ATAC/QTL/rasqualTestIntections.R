library("devtools")
library("plyr")
library("dplyr")
library("ggplot2")
library("purrr")
load_all("../seqUtils/")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")

#Import data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
gene_name_map = dplyr::select(atac_list$gene_metadata, gene_id, gene_name)

#Import minimal p-values
rasqual_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Rasqual selected p-values
rasqual_selected_results = readRDS("results/ATAC/QTLs/rasqual_selected_pvalues.rds")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Calculate R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, filtered_vcf$genotypes, .8)

#### Test for interactions ####
#covariate_names = c("sex_binary","short_long_ratio","assigned_frac","mt_frac")
covariate_names = c("sex_binary", "cqn_PC1", "cqn_PC2", "cqn_PC3")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Run model
interaction_results = testMultipleInteractions(filtered_pairs, trait_matrix = atac_list$cqn, 
                    sample_metadata = atac_list$sample_metadata, filtered_vcf, formula_qtl, formula_interaction, id_field_separator = "-")
interaction_df = postProcessInteractionPvalues(interaction_results, id_field_separator = "-")
saveRDS(interaction_df, "results/ATAC/QTLs/rasqual_interaction_results.rds")
interaction_df = readRDS("results/ATAC/QTLs/rasqual_interaction_results.rds")
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)


#Make a Q-Q plot for the interaction p-values
qq_df = dplyr::mutate(interaction_df, p_eigen = p_nominal) %>% addExpectedPvalue()
qq_plot = ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_nominal,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")
ggsave("figures/supplementary/caQTL_interaction_Q-Q_plot_tpm.pdf", plot = qq_plot, width = 4, height = 4)

#Quantify genomic inflation
lambda = median(qchisq(1-interaction_df$p_nominal,1))/qchisq(0.5,1)

#### Cluster effect sizes ####
#Extract effect sizes for all gene-snp pairs from RASQUAL data
beta_list = extractAndProcessBetas(dplyr::select(interaction_hits, gene_id, snp_id), rasqual_selected_results, "naive")

#Find QTLs that appear
set.seed(42)
appear_qtls = dplyr::filter(beta_list$beta_summaries, abs(naive) <= 0.59, max_abs_diff >= 0.59, max_abs_beta >= 0.59)
appear_betas = dplyr::semi_join(beta_list$beta_summaries[,1:6], appear_qtls, by = c("gene_id", "snp_id"))
appear_clusters = clusterBetasKmeans(appear_betas, 6) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_list$beta_df, by = c("gene_id", "snp_id"))

#Calculate mean effect size in each cluster and condition
appear_cluster_means = calculateClusterMeans(appear_clusters)
cluster_sizes = calculateClusterSizes(appear_clusters)

#Reorder clusters
cluster_reorder = data_frame(cluster_id = c(1,2,3,4,5,6), new_cluster_id = c(5,1,2,3,6,4))

#Make heatmap of effect sizes
appear_betas = appear_clusters %>% dplyr::group_by(gene_id, snp_id) %>% 
  dplyr::mutate(beta_scaled = beta/max(beta)) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>% 
  dplyr::semi_join(cluster_sizes, by = "cluster_id") %>% 
  ungroup() %>%
  dplyr::left_join(cluster_reorder, by = "cluster_id") %>%
  dplyr::left_join(figureNames(), by = "condition_name") #Add figure names


#Make a plot of effect sizes
ylabel = paste(sum(cluster_sizes$count), "caQTLs")
effect_size_heatmap = ggplot(appear_betas, aes(x = figure_name, y = gene_name, fill = beta_scaled)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ylab(ylabel) + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Relative effect", midpoint = 0) +
  theme_light() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.x = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))
ggsave("figures/main_figures/caQTLs_appear_kmeans_heatmap.png",effect_size_heatmap, width = 3.5, height = 4)


#Find QTLs that disappear
#Look for QTLs that disappear after stimulation
disappear_qtls = dplyr::filter(beta_list$beta_summaries, abs(naive) > 0.59, max_abs_diff >= 0.59, min_abs_beta <= 0.59)
disappear_betas = dplyr::semi_join(beta_list$beta_summaries[,1:6], disappear_qtls, by = c("gene_id", "snp_id"))
disappear_clusters = clusterBetasKmeans(disappear_betas, 7) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_list$beta_df, by = c("gene_id", "snp_id"))
disappear_plot = ggplot(disappear_clusters, aes(x = condition_name, y = abs(beta), group = paste(gene_id, snp_id))) + 
  geom_line() + facet_wrap(~cluster_id)

#Calculate mean effect size in each cluster and condition
disappear_cluster_sizes = calculateClusterSizes(disappear_clusters, selected_condition = "IFNg_SL1344") %>% 
  dplyr::filter(count > 50) #Remove small clusters
disappear_cluster_means = calculateClusterMeans(disappear_clusters) %>%
  dplyr::semi_join(disappear_cluster_sizes, by = "cluster_id")

disappear_betas = disappear_clusters %>% dplyr::group_by(gene_id, snp_id) %>% dplyr::mutate(beta_scaled = beta/max(beta)) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>%
  dplyr::semi_join(disappear_cluster_sizes, by = "cluster_id") %>% 
  ungroup() %>%
  dplyr::left_join(figureNames(), by = "condition_name")

ylabel = paste(sum(disappear_cluster_sizes$count), "caQTLs")
dis_effect_size_heatmap = ggplot(disappear_betas, aes(x = figure_name, y = gene_name, fill = beta_scaled)) + 
  facet_grid(cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ylab(ylabel) + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Relative effect", midpoint = 0) +
  theme_grey() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), 
        axis.title.x = element_blank()) +
  theme(panel.spacing = unit(0.2, "lines"))
ggsave("figures/supplementary/caQTLs_disappear_kmeans_heatmap.png",dis_effect_size_heatmap, width = 3.5, height = 4)

#Export variable qtls
variable_qtls = list(appear = appear_betas, disappear = disappear_betas)
saveRDS(variable_qtls, "results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")
