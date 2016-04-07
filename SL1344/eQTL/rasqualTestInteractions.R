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

#Load p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_min_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = plyr::ldply(rasqual_min_hits, .id = "condition_name") %>% dplyr::arrange(p_nominal)
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
interaction_results = testMultipleInteractions(tbl_df(filtered_pairs), combined_expression_data$cqn, combined_expression_data$sample_metadata, filtered_vcf, formula_qtl, formula_interaction)
interaction_df = postProcessInteractionPvalues(interaction_results)
saveRDS(interaction_df, "results/SL1344/eQTLs/SL1344_interaction_pvalues.rds")
interaction_df = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues.rds")
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

#Make heatmap of effect sizes
appear_betas = appear_clusters %>% dplyr::group_by(gene_id, snp_id) %>% dplyr::mutate(beta_scaled = beta/max(beta)) %>%
  dplyr::left_join(gene_name_map, by = "gene_id")
effect_size_heatmap = ggplot(appear_betas, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + 
  facet_grid(cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Beta", midpoint = 0) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave("results/SL1344/eQTLs/properties/eQTLs_appear_kmeans_heatmap.pdf",effect_size_heatmap, width = 5, height = 7)


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

disappear_betas = disappear_clusters %>% dplyr::group_by(gene_id, snp_id) %>% dplyr::mutate(beta_scaled = beta/max(beta)) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>%
  dplyr::semi_join(disappear_cluster_sizes, by = "cluster_id")
dis_effect_size_heatmap = ggplot(disappear_betas, aes(x = condition_name, y = gene_name, fill = beta_scaled)) + 
  facet_grid(cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Beta", midpoint = 0) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave("results/SL1344/eQTLs/properties/eQTLs_disappear_kmeans_heatmap.pdf",dis_effect_size_heatmap, width = 5, height = 5)


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


#### GWAS overlaps ####
#Import GWAS catalog
filtered_catalog = readRDS("annotations/gwas_catalog_v1.0.1-downloaded_2016-03-02.filtered.rds")

#GWAS overlaps for QTLs that appear
qtl_table = dplyr::select(interaction_hits, gene_id, snp_id)
interaction_olaps = findGWASOverlaps(qtl_table, filtered_catalog, vcf_file)
interaction_gwas_hits = dplyr::left_join(interaction_olaps, gene_name_map, by = "gene_id") %>%
  dplyr::select(gene_name, snp_id, gwas_snp_id, R2, trait, gwas_pvalue)
write.table(interaction_gwas_hits, "results/SL1344/eQTLs/appear_gwas_overlaps.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#All GWAS overlaps
all_olaps = findGWASOverlaps(filtered_pairs, filtered_catalog, vcf_file, min_r2 = 0.6)
all_gwas_hits = dplyr::left_join(all_olaps, gene_name_map, by = "gene_id") %>%
  dplyr::select(gene_name, gene_id, snp_id, gwas_snp_id, R2, trait, gwas_pvalue)
write.table(all_gwas_hits, "results/SL1344/eQTLs/all_gwas_overlaps.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Rank traits by overlap size
ranked_traits = rankTraitsByOverlapSize(dplyr::filter(all_gwas_hits, R2 > 0.78), filtered_catalog, min_overlap = 4)
write.table(ranked_traits, "results/SL1344/eQTLs/relative_gwas_overlaps_R08.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Compare ATAC and RNA overlaps
atac_olaps = readr::read_delim("../macrophage-chromatin/results/ATAC/QTLs/all_gwas_overlaps_R08.txt", delim = "\t")
rna_overlaps = dplyr::filter(all_gwas_hits, R2 > 0.8)

#Look at the overlap between GWAS overlaps
atac_unique_snps = unique(atac_olaps$gwas_snp_id)
rna_unique_snps = unique(rna_overlaps$gwas_snp_id)
intersect(atac_unique_snps, rna_unique_snps)





#Explore multiple SNPs per gene
#Fetch p-values from rasqual output
gene_df = data_frame(gene_id = "ENSG00000091490")
gene_ranges = constructGeneRanges(gene_df, combined_expression_data$gene_metadata, cis_window = 5e5)
tabix_data = tabixFetchGenes(gene_ranges, "databases/SL1344/IFNg_500kb.sorted.txt.gz")[[1]] %>%
  dplyr::arrange(p_nominal) %>%
  addAssociationPosterior(86) %>%
  dplyr::filter(lABF > max(lABF)/3)
unique_snps = filterGeneR2(tabix_data, vcf_file$genotypes, 0.2)
plot(tabix_data$pos, tabix_data$lABF)



tabix_files = list(naive = "databases/SL1344/naive_500kb.sorted.txt.gz", IFNg = "databases/SL1344/IFNg_500kb.sorted.txt.gz",
                   SL1344 = "databases/SL1344/SL1344_500kb.sorted.txt.gz", IFNg_SL1344 = "databases/SL1344/IFNg_SL1344_500kb.sorted.txt.gz")
coords = lapply(tabix_files, function(x, gene_ranges) {
  tabixFetchGenes(gene_ranges, x)[[1]]
}, gene_ranges)
beta_list1 = extractAndProcessBetas(dplyr::select(unique_snps, gene_id, snp_id), coords, "naive")
beta_list2 = extractAndProcessBetas(dplyr::select(head(tabix_data,4), gene_id, snp_id), coords, "naive")

#PTK2B
gene_df = data_frame(gene_id = "ENSG00000120899")
gene_ranges = constructGeneRanges(gene_df, combined_expression_data$gene_metadata, cis_window = 5e5)
tabix_data = tabixFetchGenes(gene_ranges, "databases/SL1344/IFNg_SL1344_500kb.sorted.txt.gz")[[1]] %>%
  dplyr::arrange(p_nominal) %>%
  addAssociationPosterior(86) %>%
  dplyr::filter(lABF > max(lABF)/3)

unique_snps = filterGeneR2(tabix_data, vcf_file$genotypes, 0.2)



tabix_files = list(naive = "databases/SL1344/naive_500kb.sorted.txt.gz", IFNg = "databases/SL1344/IFNg_500kb.sorted.txt.gz",
                   SL1344 = "databases/SL1344/SL1344_500kb.sorted.txt.gz", IFNg_SL1344 = "databases/SL1344/IFNg_SL1344_500kb.sorted.txt.gz")
coords = lapply(tabix_files, function(x, gene_ranges) {
  tabixFetchGenes(gene_ranges, x)[[1]]
}, gene_ranges)
beta_list1 = extractAndProcessBetas(dplyr::select(unique_snps, gene_id, snp_id), coords, "naive")
beta_list2 = extractAndProcessBetas(dplyr::select(head(tabix_data,4), gene_id, snp_id), coords, "naive")



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
