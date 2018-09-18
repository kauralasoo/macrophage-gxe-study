library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

#Functions
quantileNormaliseBeta <- function(beta){
  m = mean(beta)
  new_beta = quantileNormaliseVector(beta - m) + m
}

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

#Import permutation p-values for the linear model
empirical_pvalues_lm = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues_lm.empirical.rds")
empirical_hits_lm = dplyr::filter(empirical_pvalues_lm, p_fdr < 0.1)

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
write.table(interaction_df, "figures/tables/RNA_eQTL_lme4_interaction_test.txt", sep = "\t", quote = FALSE, row.names = FALSE)
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)

#Import permutation p-values for the linear model
empirical_pvalues_lme4 = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues_lme4.empirical.rds")
empirical_hits_lme4 = dplyr::filter(empirical_pvalues_lme4, p_fdr < 0.1)


#Make a QQ-plot
qq_df = dplyr::mutate(interaction_df, p_eigen = p_nominal) %>% addExpectedPvalue()
qq_plot = ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_nominal,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")

#Comapre lm and lme4 p-values
lme4_pval = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues_lme4.rds") %>%
  dplyr::transmute(gene_id, snp_id, p_lme4 = p_nominal)
lm_pval = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues.rds") %>%
  dplyr::transmute(gene_id, snp_id, p_lm = p_nominal)
joint_p = dplyr::left_join(lme4_pval, lm_pval, by = c("gene_id", "snp_id"))
lme4_plot = ggplot(joint_p, aes(x = -log(p_lm, 10), y = -log(p_lme4, 10))) + geom_point() +
  theme_light() +
  xlab("linear model p-value") +
  ylab("linear mixed model p-value")
ggsave("figures/supplementary/eQTL_lm_vs_lme4_pvalues.png",lme4_plot, width = 5, height = 5)


#Extract effect sizes for all gene-snp pairs from RASQUAL data
beta_list = extractAndProcessBetas(dplyr::select(interaction_hits, gene_id, snp_id), rasqual_selected_pvalues, "naive")
beta_list$beta_summaries = dplyr::mutate(beta_list$beta_summaries, max_naive_ratio = max_abs_beta/abs(naive))

#Find QTLs that appear
set.seed(42)
#appear_qtls = dplyr::filter(beta_list$beta_summaries, abs(naive) <= 0.59, max_abs_beta - abs(naive) >= 0.32)
appear_qtls = dplyr::filter(beta_list$beta_summaries, (abs(naive) <= 0.59 & max_abs_diff >= 0.59 & max_abs_beta >= 0.59) | 
                              gene_id == "ENSG00000144227")
appear_betas = dplyr::semi_join(beta_list$beta_summaries[,1:6], appear_qtls, by = c("gene_id", "snp_id"))
appear_clusters = clusterBetasKmeans(appear_betas, 6) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_list$beta_df, by = c("gene_id", "snp_id"))

#Keep only those genes that survive permutation testing
appear_clusters = dplyr::semi_join(appear_clusters, empirical_hits_lme4, by = c("gene_id", "snp_id"))

#Reorder clusters
cluster_reorder = data_frame(cluster_id = c(1,2,3,4,5,6), new_cluster_id = c(6,3,4,1,2,5))

#Make heatmap of effect sizes
appear_betas = appear_clusters %>% 
  dplyr::group_by(gene_id, snp_id) %>% 
  dplyr::mutate(beta_scaled = beta/max(beta)) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>%
  dplyr::left_join(cluster_reorder, by = "cluster_id") %>% #Reorder clusters
  dplyr::left_join(figureNames(), by = "condition_name") %>% #Add figure names %>%
  dplyr::group_by(gene_id, snp_id) %>% 
  dplyr::arrange(gene_id, snp_id, -abs(beta)) %>% 
  dplyr::mutate(max_condition = condition_name[1]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_quantile = quantileNormaliseBeta(beta))

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
  theme(legend.title = element_text(angle = 90)) + 
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.x = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(strip.text.y = element_text(colour = "grey10"), strip.background = element_rect(fill = "grey85"))
ggsave("figures/main_figures/eQTLs_appear_kmeans_heatmap.png",effect_size_heatmap, width = 3, height = 4)


#Count conditions with maximal effect sizes
max_condition_table = dplyr::select(appear_betas, gene_id, snp_id, max_condition) %>% 
  unique() %>%
  dplyr::group_by(max_condition) %>%
  dplyr::summarise(max_count = length(max_condition)) %>%
  dplyr::rename(condition_name = max_condition) %>%
  dplyr::left_join(figureNames())

max_effect_plot = ggplot(max_condition_table, aes(x = figure_name, y = max_count)) + 
  geom_bar(stat = "identity") + 
  xlab("Condition") + 
  theme_light() +
  ylab("Response eQTL count")
ggsave("figures/supplementary/eQTL_max_condition_count.pdf", plot = max_effect_plot, width = 3, height = 3)

#Make a quantile normalised plot
effect_size_heatmap = ggplot(appear_betas, aes(x = figure_name, y = gene_name, fill = beta_quantile)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ylab(ylabel) + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Relative effect", midpoint = 0) +
  theme_light() +
  theme(legend.title = element_text(angle = 90)) + 
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.x = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(strip.text.y = element_text(colour = "grey10"), strip.background = element_rect(fill = "grey85"))
ggsave("figures/main_figures/eQTLs_appear_kmeans_heatmap_quantile.png",effect_size_heatmap, width = 3, height = 4)

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
disappear_qtls = dplyr::filter(beta_list$beta_summaries, (abs(naive) >= 0.59 & max_abs_diff >= 0.59 & min_abs_beta <= 0.59))
disappear_betas = dplyr::semi_join(beta_list$beta_summaries[,1:6], disappear_qtls, by = c("gene_id", "snp_id"))
disappear_clusters = clusterBetasKmeans(disappear_betas, 7) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_list$beta_df, by = c("gene_id", "snp_id"))
disappear_clusters = dplyr::semi_join(disappear_clusters, empirical_hits_lme4, by = c("gene_id", "snp_id"))


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
  dplyr::left_join(figureNames(), by = "condition_name") %>% #Add figure names
  dplyr::ungroup()

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


disappear_means_plot = ggplot(disappear_cluster_means, aes(x = figure_name, y = beta_mean, group = cluster_id)) + 
  geom_point() + geom_line() + facet_wrap(~cluster_id) + 
  geom_errorbar(aes(ymin = beta_mean - beta_sd, ymax = beta_mean + beta_sd, width = .2)) + 
  geom_text(aes(label = paste("n = ", count, sep = ""), x = condition_name, y = 1.5), data = disappear_cluster_sizes) +
  xlab("Condition") + 
  ylab("Log2 fold-change")

#Export clustering results
variable_qtls = list(appear = appear_betas, disappear = disappear_betas)
saveRDS(variable_qtls, "results/SL1344/eQTLs/appear_disappear_eQTLs.rds")

#Export gene lists to disk
variable_qtls = readRDS("results/SL1344/eQTLs/appear_disappear_eQTLs.rds")
write.table(variable_qtls$appear, "figures/tables/Figure_2a_response_eQTLs.txt", sep = "\t", quote = FALSE, row.names = FALSE)






#Compare Wald test to LRT
naive_ifng_data = extractConditionFromExpressionList(c("naive","IFNg"), combined_expression_data)
naive_sl1344_data = extractConditionFromExpressionList(c("naive","SL1344"), combined_expression_data)
naive_ifng_sl1344_data = extractConditionFromExpressionList(c("naive","IFNg_SL1344"), combined_expression_data)

#Perform interaction test using Wald test
interaction_results = testMultipleInteractions(tbl_df(filtered_pairs), naive_ifng_data$cqn, 
                                               naive_ifng_data$sample_metadata, 
                                               filtered_vcf, formula_qtl, formula_interaction, id_field_separator = "-", return_value = "model")
wald_res = purrr::map_df(interaction_results, ~data_frame(p_wald = (summary(.$interaction_model) %>% coef())[11,4]), .id = "id") %>%
  tidyr::separate(id, into = c("gene_id", "snp_id"), sep = "-") 

interaction_results = testMultipleInteractions(tbl_df(filtered_pairs), naive_ifng_data$cqn, 
                                               naive_ifng_data$sample_metadata, 
                                               filtered_vcf, formula_qtl, formula_interaction, id_field_separator = "-")
anova_res = postProcessInteractionPvalues(interaction_results, id_field_separator = "-") %>%
  dplyr::transmute(gene_id, snp_id, p_lrt = p_nominal)


anova_lm_plot = ggplot(df, aes(x = -log(p_lrt,10), y = -log(df$p_wald,10))) + 
  geom_point() +
  xlab("Anova p-value") +
  ylab("summary(lm) p-value") + 
  theme_light()
ggsave("figures/supplementary/anova_vs_lm_pvalues.png", plot = anova_lm_plot, width = 4, height = 4)

