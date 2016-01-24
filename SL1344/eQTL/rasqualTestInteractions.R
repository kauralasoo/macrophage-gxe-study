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

rasqual_selected_pvalues = readRDS("results/SL1344/eQTLs/rasqual_selected_pvalues.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Calculate R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = vcf_file$snpspos
snps_pos2 = dplyr::filter(snps_pos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos2, genotypes = genotypes)

#Naive vs IFNg
condition_pair = c("naive", "IFNg")
covariate_names = c("sex_binary", "ng_ul_mean","macrophage_diff_days","rna_auto", "max_purity_filtered", "harvest_stimulation_days",
                    "PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6")


#Calculate all pairwise interactions
ifng_interactions = testInterctionsBetweenPairs(c("naive", "IFNg"), rasqual_min_hits, combined_expression_data, covariate_names, filtered_vcf) %>%
  dplyr::filter(abs_beta_min <= 0.59, abs(beta_diff) >= 0.59)
sl1344_interactions = testInterctionsBetweenPairs(c("naive", "SL1344"), rasqual_min_hits, combined_expression_data, covariate_names,filtered_vcf) %>%
  dplyr::filter(abs_beta_min <= 0.59, abs(beta_diff) >= 0.59)
ifng_sl1344_interactions = testInterctionsBetweenPairs(c("naive", "IFNg_SL1344"), rasqual_min_hits, combined_expression_data, covariate_names, filtered_vcf) %>%
  dplyr::filter(abs_beta_min <= 0.59, abs(beta_diff) >= 0.59)

#Make plots
makeMultiplePlots(ifng_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_IFNg/", 7,7)
makeMultiplePlots(sl1344_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_SL1344/", 7,7)
makeMultiplePlots(ifng_sl1344_interactions, combined_expression_data$cqn,filtered_vcf$genotypes,combined_expression_data$sample_metadata,combined_expression_data$gene_metadata) %>%
  savePlots("results/SL1344/eQTLs/interaction_plots/naive_vs_IFNg_SL1344/", 7,7)

