library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")


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

#Naive vs IFNg
covariate_names = c("sex_binary", "ng_ul_mean","macrophage_diff_days","rna_auto", "max_purity_filtered", "harvest_stimulation_days",
                    "PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6")


#Calculate all pairwise interactions
ifng_interactions = testInterctionsBetweenPairs(c("naive", "IFNg"), rasqual_min_hits, combined_expression_data, covariate_names, filtered_vcf) %>%
  dplyr::filter(abs_beta_min <= 0.59, abs(beta_diff) >= 0.59)
sl1344_interactions = testInterctionsBetweenPairs(c("naive", "SL1344"), rasqual_min_hits, combined_expression_data, covariate_names,filtered_vcf) %>%
  dplyr::filter(abs_beta_min <= 0.59, abs(beta_diff) >= 0.59)
ifng_sl1344_interactions = testInterctionsBetweenPairs(c("naive", "IFNg_SL1344"), rasqual_min_hits, combined_expression_data, covariate_names, filtered_vcf) %>%
  dplyr::filter(abs_beta_min <= 0.59, abs(beta_diff) >= 0.59)
all_interactions = rbind(ifng_sl1344_interactions, sl1344_interactions, ifng_interactions) %>%
  dplyr::select(gene_id, snp_id) %>% unique()

#Extract effect sizes for all gene-snp pairs from RASQUAL data
beta_matrix = extractBetasFromList(all_interactions, rasqual_selected_pvalues) %>% ungroup()

#Calculate maximum absolute diff in betas between conditions
beta_diff_matrix = dplyr::select(beta_matrix, -gene_id, -snp_id)
beta_diff_matrix = beta_diff_matrix - beta_diff_matrix$naive
beta_diff_matrix$max_abs_diff = apply(abs(beta_diff_matrix), 1, max)
beta_diff_matrix = cbind(beta_matrix[,c("gene_id", "snp_id")], beta_diff_matrix)
  
#Convert Beta matrix into a df
beta_df = tidyr::gather(beta_matrix, condition, beta, naive:IFNg_SL1344)
beta_summaries = dplyr::group_by(beta_df, gene_id, snp_id) %>% 
  dplyr::summarise(max_abs_beta = max(abs(beta)), min_abs_beta = min(abs(beta))) %>%
  dplyr::left_join(dplyr::select(beta_diff_matrix, gene_id, snp_id, max_abs_diff), by = c("gene_id","snp_id"))
beta_summaries_matrix = dplyr::left_join(beta_matrix, beta_summaries)

#Find QTLs that appear
appear_qtls = dplyr::filter(beta_summaries_matrix, abs(naive) <= 0.59, max_abs_diff >= 0.59, max_abs_beta >= 0.59)
appear_betas = dplyr::semi_join(beta_matrix, appear_qtls, by = c("gene_id", "snp_id"))
appear_clusters = clusterBetasKmeans(appear_betas, 5) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_df, by = c("gene_id", "snp_id"))
appear_plot = ggplot(appear_clusters, aes(x = condition, y = abs(beta), group = paste(gene_id, snp_id))) + geom_line() + facet_wrap(~cluster_id)
ggsave("results/SL1344/eQTLs/properties/eQTLs_appear_kmeans.pdf",appear_plot, width = 10, height = 10)

#Find QTLs that disappear
#Look for QTLs that disappear after stimulation
disappear_qtls = dplyr::filter(beta_summaries_matrix, abs(naive) > 0.59, max_abs_diff >= 0.59, min_abs_beta <= 0.59)
disappear_betas = dplyr::semi_join(beta_matrix, disappear_qtls, by = c("gene_id", "snp_id"))
disappear_clusters = clusterBetasKmeans(disappear_betas, 7) %>% dplyr::select(gene_id, snp_id, cluster_id) %>%
  dplyr::left_join(beta_df, by = c("gene_id", "snp_id"))
disappear_plot = ggplot(disappear_clusters, aes(x = condition, y = abs(beta), group = paste(gene_id, snp_id))) + geom_line() + facet_wrap(~cluster_id)
ggsave("results/SL1344/eQTLs/properties/eQTLs_disappear_kmeans.pdf",disappear_plot, width = 10, height = 10)


#Make plots
makeMultiplePlots(ifng_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_IFNg/", 7,7)
makeMultiplePlots(sl1344_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_SL1344/", 7,7)
makeMultiplePlots(ifng_sl1344_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_IFNg_SL1344/", 7,7)

