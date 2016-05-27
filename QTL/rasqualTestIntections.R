library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")
library("ggplot2")
library("purrr")

#Import data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
gene_name_map = dplyr::select(atac_list$gene_metadata, gene_id, gene_name)

#Import minimal p-values
min_pvalue_list = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(min_pvalue_list, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::group_by(gene_id) %>% dplyr::arrange(p_nominal) %>% ungroup()
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Rasqual selected p-values
rasqual_selected_results = readRDS("results/ATAC/QTLs/rasqual_selected_results.rds")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Calculate R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, filtered_vcf$genotypes, .8)

#### Test for interactions ####
covariate_names = c("sex_binary", "cqn_PC1", "cqn_PC2", "cqn_PC3")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Run model
interaction_results = testMultipleInteractions(filtered_pairs, trait_matrix = atac_list$cqn, 
                    sample_metadata = atac_list$sample_metadata, filtered_vcf, formula_qtl, formula_interaction)
interaction_df = postProcessInteractionPvalues(interaction_results)
saveRDS(interaction_df, "results/ATAC/QTLs/rasqual_interaction_results.rds")
interaction_df = readRDS("results/ATAC/QTLs/rasqual_interaction_results.rds")
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)


#### Cluster effect sizes ####
#Extract effect sizes for all gene-snp pairs from RASQUAL data
beta_list = extractAndProcessBetas(dplyr::select(interaction_hits, gene_id, snp_id), rasqual_selected_results, "naive")

#Find QTLs that appear
set.seed(42)
appear_qtls = dplyr::filter(beta_list$beta_summaries, abs(naive) <= 0.59, max_abs_diff >= 0.59, max_abs_beta >= 0.59)
appear_betas = dplyr::semi_join(beta_list$beta_summaries[,1:6], appear_qtls, by = c("gene_id", "snp_id"))
appear_clusters = clusterBetasKmeans(appear_betas, 6) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_list$beta_df, by = c("gene_id", "snp_id"))
appear_plot = ggplot(appear_clusters, aes(x = condition_name, y = beta, group = paste(gene_id, snp_id))) + 
  geom_line() + facet_wrap(~cluster_id)

#Calculate mean effect size in each cluster and condition
appear_cluster_means = calculateClusterMeans(appear_clusters)
cluster_sizes = calculateClusterSizes(appear_clusters)

#Make heatmap of effect sizes
appear_betas = appear_clusters %>% dplyr::group_by(gene_id, snp_id) %>% dplyr::mutate(beta_scaled = beta/max(beta)) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>% 
  dplyr::semi_join(cluster_sizes, by = "cluster_id") %>% 
  ungroup()

#Reorder clusters
cluster_reorder = data_frame(cluster_id = c(1,2,3,4,5,6), new_cluster_id = c(1,6,4,3,5,2))
appear_betas = dplyr::left_join(appear_betas, cluster_reorder, by = "cluster_id")

#Make a plot of effect sizes
ylabel = paste(sum(cluster_sizes$count), "caQTLs")
effect_size_heatmap = ggplot(appear_betas, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ylab(ylabel) + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Relative effect", midpoint = 0) +
  theme_grey() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), 
        axis.title.x = element_blank(), axis.text.x=element_text(angle = 15)) +
  theme(panel.margin = unit(0.2, "lines"))
ggsave("results/ATAC/QTLs/properties/caQTLs_appear_kmeans_heatmap.png",effect_size_heatmap, width = 4.5, height = 5.5)


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
  dplyr::semi_join(disappear_cluster_sizes, by = "cluster_id") %>% ungroup()

ylabel = paste(sum(disappear_cluster_sizes$count), "caQTLs")
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
ggsave("results/ATAC/QTLs/properties/caQTLs_disappear_kmeans_heatmap.png",dis_effect_size_heatmap, width = 4.5, height = 4)

#Export variable qtls
variable_qtls = list(appear = appear_betas, disappear = disappear_betas)
saveRDS(variable_qtls, "results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")

#Plot example
plotEQTL("ATAC_peak_26312","rs789642", atac_list$cqn, vcf_file$genotypes, atac_list$sample_metadata, atac_list$gene_metadata)


#### GWAS overlaps ####
#Import GWAS catalog
filtered_catalog = readRDS("../macrophage-gxe-study/annotations/gwas_catalog_v1.0.1-downloaded_2016-03-02.filtered.rds")

#GWAS overlaps for QTLs that appear
qtl_table = dplyr::select(interaction_hits, gene_id, snp_id)
interaction_olaps = findGWASOverlaps(qtl_table, filtered_catalog, vcf_file)
interaction_gwas_hits = dplyr::left_join(interaction_olaps, gene_name_map, by = "gene_id") %>%
  dplyr::select(gene_name, snp_id, gwas_snp_id, R2, trait, gwas_pvalue)
write.table(interaction_gwas_hits, "results/ATAC/eQTLs/appear_gwas_overlaps.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#All GWAS overlaps
all_olaps = findGWASOverlaps(filtered_pairs, filtered_catalog, vcf_file, min_r2 = 0.8)
all_gwas_hits = dplyr::left_join(all_olaps, gene_name_map, by = "gene_id") %>%
  dplyr::select(gene_name, gene_id, snp_id, gwas_snp_id, R2, trait, gwas_pvalue)
write.table(all_gwas_hits, "results/ATAC/QTLs/all_gwas_overlaps_R08.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Rank traits by overlap size
ranked_traits = rankTraitsByOverlapSize(dplyr::filter(all_gwas_hits, R2 > 0.8), filtered_catalog, min_overlap = 4)
write.table(ranked_traits, "results/ATAC/QTLs/relative_gwas_overlaps_R08.txt", row.names = FALSE, sep = "\t", quote = FALSE)


