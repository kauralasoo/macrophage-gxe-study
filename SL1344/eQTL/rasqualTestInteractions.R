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
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Import summary stats for each pair in each condition
rasqual_selected_pvalues = readRDS("results/SL1344/eQTLs/rasqual_selected_pvalues.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Calculate R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

#Naive vs IFNg
covariate_names = c("PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6", "sex_binary")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                      paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(tbl_df(filtered_pairs), combined_expression_data$cqn, 
                                               combined_expression_data$sample_metadata, 
                                               filtered_vcf, formula_qtl, formula_interaction, id_field_separator = "-")
interaction_df = postProcessInteractionPvalues(interaction_results, id_field_separator = "-")
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

#Use a paired design to test for interaction
covariate_names = c("PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6", "sex_binary")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(tbl_df(filtered_pairs), combined_expression_data$cqn, 
                                               combined_expression_data$sample_metadata, 
                                               filtered_vcf, formula_qtl, formula_interaction, id_field_separator = "-", lme4 = TRUE)
interaction_df = postProcessInteractionPvalues(interaction_results, id_field_separator = "-")
saveRDS(interaction_df, "results/SL1344/eQTLs/SL1344_interaction_pvalues_lme4.rds")
interaction_df = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues_lme4.rds")
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)

qq_df = dplyr::mutate(interaction_df, p_eigen = p_nominal) %>% addExpectedPvalue()
qq_plot = ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_nominal,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")




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
cluster_reorder = data_frame(cluster_id = c(1,2,3,4,5,6), new_cluster_id = c(1,6,5,4,3,2))

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
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(strip.text.y = element_text(colour = "grey10"), strip.background = element_rect(fill = "grey85"))
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

#Find QTLs that disappear
#Look for QTLs that disappear after stimulation
set.seed(41)
disappear_qtls = dplyr::filter(beta_list$beta_summaries, abs(naive) > 0.59, abs(naive) - min_abs_beta > 0.32)
disappear_betas = dplyr::semi_join(beta_list$beta_summaries[,1:6], disappear_qtls, by = c("gene_id", "snp_id"))
disappear_clusters = clusterBetasKmeans(disappear_betas, 7) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_list$beta_df, by = c("gene_id", "snp_id"))
disappear_plot = ggplot(disappear_clusters, aes(x = condition_name, y = abs(beta), group = paste(gene_id, snp_id))) + 
  geom_line() + facet_wrap(~cluster_id)


#Calculate mean effect size in each cluster and condition
disappear_cluster_sizes = calculateClusterSizes(disappear_clusters, selected_condition = "IFNg_SL1344") %>% 
  dplyr::filter(count > 10)
disappear_cluster_means = calculateClusterMeans(disappear_clusters) %>%
  dplyr::semi_join(disappear_cluster_sizes, by = "cluster_id")

disappear_betas = disappear_clusters %>% 
  dplyr::group_by(gene_id, snp_id) %>% 
  dplyr::mutate(beta_scaled = beta/max(beta)) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>%
  dplyr::semi_join(disappear_cluster_sizes, by = "cluster_id") %>%
  dplyr::left_join(figureNames(), by = "condition_name") #Add figure names

disappear_count = disappear_cluster_sizes$count %>% sum()
ylabel = paste(disappear_count, "eQTLs")
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
ggsave("figures/supplementary/eQTLs_disappear_kmeans_heatmap.png",dis_effect_size_heatmap, width = 3.5, height = 3)


disappear_means_plot = ggplot(disappear_cluster_means, aes(x = condition_name, y = beta_mean, group = cluster_id)) + 
  geom_point() + geom_line() + facet_wrap(~cluster_id) + 
  geom_errorbar(aes(ymin = beta_mean - beta_sd, ymax = beta_mean + beta_sd, width = .2)) + 
  geom_text(aes(label = paste("n = ", count, sep = ""), x = condition_name, y = 1.5), data = disappear_cluster_sizes) +
  xlab("Condition") + 
  ylab("Log2 fold-change")

#Export clustering results
variable_qtls = list(appear = appear_betas, disappear = disappear_betas)
saveRDS(variable_qtls, "results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")


