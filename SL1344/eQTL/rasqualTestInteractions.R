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
combined_expression_data$gene_metadata = dplyr::rename(combined_expression_data$gene_metadata, 
                                                       chr = chromosome_name, start = start_position, end = end_position)
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                  levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))

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
beta_matrix = extractBetasFromList(dplyr::select(interaction_hits, gene_id, snp_id), rasqual_selected_pvalues) %>% ungroup()

#Calculate maximum absolute diff in betas between conditions
beta_diff_matrix = dplyr::select(beta_matrix, -gene_id, -snp_id)
beta_diff_matrix = beta_diff_matrix - beta_diff_matrix$naive
beta_diff_matrix$max_abs_diff = apply(abs(beta_diff_matrix), 1, max)
beta_diff_matrix = cbind(beta_matrix[,c("gene_id", "snp_id")], beta_diff_matrix)
  
#Convert Beta matrix into a df and revert signs to the sign of the maximal effect size
beta_df = tidyr::gather(beta_matrix, condition_name, beta, naive:IFNg_SL1344)
beta_correct_sign = betaCorrectSign(beta_df)

#Convert to matrix
beta_correct_sign_matrix = tidyr::spread(beta_correct_sign, condition_name, beta)

#Merge summary stats
beta_summaries = dplyr::group_by(beta_correct_sign, gene_id, snp_id) %>% 
  dplyr::summarise(max_abs_beta = max(abs(beta)), min_abs_beta = min(abs(beta)))%>%
  ungroup() %>%
  dplyr::left_join(dplyr::select(beta_diff_matrix, gene_id, snp_id, max_abs_diff), by = c("gene_id","snp_id"))
beta_summaries_matrix = dplyr::left_join(beta_correct_sign_matrix, beta_summaries, by = c("gene_id","snp_id"))

#Find QTLs that appear
appear_qtls = dplyr::filter(beta_summaries_matrix, abs(naive) <= 0.59, max_abs_diff >= 0.59, max_abs_beta >= 0.59)
appear_betas = dplyr::semi_join(beta_correct_sign_matrix, appear_qtls, by = c("gene_id", "snp_id"))
appear_clusters = clusterBetasKmeans(appear_betas, 6) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_correct_sign, by = c("gene_id", "snp_id"))
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
disappear_qtls = dplyr::filter(beta_summaries_matrix, abs(naive) > 0.59, max_abs_diff >= 0.59, min_abs_beta <= 0.59)
disappear_betas = dplyr::semi_join(beta_correct_sign_matrix, disappear_qtls, by = c("gene_id", "snp_id"))
disappear_clusters = clusterBetasKmeans(disappear_betas, 7) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_correct_sign, by = c("gene_id", "snp_id"))
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

#Extract most associated peaks for each of the eQTL genes

gene_name_map =  dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)
appear_hits = dplyr::left_join(appear_qtls, gene_name_map, by = "gene_id")

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")

#Load ATAC QTL hits
#Extract all SNP-gene pairs
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Construct GRanges object of SNP positions
selected_snps = dplyr::semi_join(filtered_vcf$snpspos, appear_hits, by = c("snpid" = "snp_id")) %>%
  dplyr::transmute(snp_id = snpid, seqnames = chr, start = pos, end = pos, strand = "*") %>%
  dataFrameToGRanges()

#Make database connections
naive_tabix = "~/projects/macrophage-chromatin/naive_100kb.sorted.unskipped.txt.gz"
IFNg_tabix = "~/projects/macrophage-chromatin/IFNg_100kb.sorted.unskipped.txt.gz"

#Fetch p-values
naive_hits = tabixFetchSNPs(selected_snps, naive_tabix)
IFNg_hits = tabixFetchSNPs(selected_snps, IFNg_tabix)

ifng_qtls = dplyr::filter(appear_clusters, cluster_id %in% c(1,3,5)) %>% 
  dplyr::select(gene_id, snp_id) %>% unique() %>%
  dplyr::left_join(beta_summaries_matrix, by = c("gene_id", "snp_id")) %>%
  dplyr::mutate(diff = IFNg - naive) %>% 
  dplyr::filter(diff >= 0.59, IFNg >= 0.59) %>%
  dplyr::group_by(gene_id) %>% dplyr::arrange(-max_abs_beta) %>%
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()

#Find candidate ATAC peaks that overlap with condition-specific eQTLs
ifng_effect_sizes = dplyr::semi_join(IFNg_hits, ifng_qtls, by = "snp_id") %>% 
  dplyr::group_by(snp_id) %>% dplyr::arrange(desc(chisq)) %>% 
  dplyr::filter(row_number() == 1) %>% dplyr::select(gene_id, snp_id, p_nominal, beta) %>% ungroup() %>%
  dplyr::semi_join(min_pvalues_hits$IFNg, by = "gene_id") %>%
  dplyr::mutate(condition_name = "IFNg") 

naive_effect_sizes = dplyr::semi_join(naive_hits, ifng_effect_sizes, by = c("gene_id", "snp_id")) %>%
  dplyr::select(gene_id, snp_id, p_nominal, beta) %>%
  dplyr::mutate(condition_name = "naive")

a = dplyr::left_join(naive_effect_sizes, ifng_effect_sizes, by = c("gene_id", "snp_id"))
a = dplyr::left_join(naive_effect_sizes, ifng_effect_sizes, by = c("gene_id", "snp_id"))
plot(a$beta.x, a$beta.y)

#Join both together
joint_effects = rbind(naive_effect_sizes, ifng_effect_sizes) %>% 
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg"))) %>%
  betaCorrectSign()
joint_effects_matrix = tidyr::spread(joint_effects, condition_name, beta) %>%
  dplyr::mutate(diff = IFNg - naive) %>%
  dplyr::mutate(phenotype = "ATAC")

#Extract the effect sizes for RNA as well
rna_effects = dplyr::semi_join(ifng_qtls, joint_effects_matrix, by = "snp_id") %>%
  dplyr::transmute(gene_id, snp_id, naive, IFNg, diff, phenotype = "RNA")

#All effects
all_effects = rbind(joint_effects_matrix, rna_effects)

#Make plot of effect size differences
plot = ggplot(all_effects, aes(x = diff)) + geom_histogram(binwidth = 0.15) + 
  facet_wrap(~phenotype, ncol = 1) + 
  geom_vline(xintercept = 0, colour = "red") +
  xlab("Effect size difference (IFNg - naive)")
ggsave("results/SL1344/eQTLs/properties/eQTLs_vs_caQTL_diff.pdf", plot, width = 6, height = 6)

#Effect size in IFNg condition
effect_plot = ggplot(all_effects, aes(x = IFNg)) + geom_histogram(binwidth = 0.15) + 
  facet_wrap(~phenotype, ncol = 1) + 
  geom_vline(xintercept = 0, colour = "red") +
  xlab("QTL effect size (IFNg)")
ggsave("results/SL1344/eQTLs/properties/eQTLs_vs_caQTL_IFNg.pdf", effect_plot, width = 6, height = 6)

ggplot(joint_effects, aes(x = condition_name, y = beta, group = snp_id)) + geom_point() + geom_line()
ggplot(a, aes(x = condition_name, y = abs(beta))) + geom_violin()


ifng_qtls1 = dplyr::filter(appear_clusters, cluster_id %in% c(4,5)) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  dplyr::semi_join(a, by = "snp_id")
ggplot(ifng_qtls1, aes(x = condition_name, y = abs(beta), group = snp_id)) + geom_point() + geom_line()

#Make plots
makeMultiplePlots(ifng_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_IFNg/", 7,7)
makeMultiplePlots(sl1344_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_SL1344/", 7,7)
makeMultiplePlots(ifng_sl1344_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_IFNg_SL1344/", 7,7)




