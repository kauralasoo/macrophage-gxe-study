library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")

#Import data
acldl_list = readRDS("results/acLDL/acLDL_combined_expression_data_covariates.rds")
acldl_list = extractConditionFromExpressionList(c("Ctrl","AcLDL"), acldl_list)

#Load p-values from disk
rasqual_min_pvalues = readRDS("results/acLDL/eQTLs/acLDL_rasqual_min_pvalues.rds")
rasqual_min_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = ldply(rasqual_min_hits, .id = "condition_name")
joint_pairs = dplyr::select(min_pvalue_df, gene_id, snp_id) %>% unique() 

#Import the VCF file
SNPRelate::snpgdsVCF2GDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz", 
                         "genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.gds", method = "copy.num.of.ref")
vcf_file = seqUtils::gdsToMatrix("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.gds")
saveRDS(vcf_file, "genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.rds")
vcf_file = readRDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.rds")

#Calculate R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Test intetactions between 
filtered_pairs = filterHitsR2(joint_pairs, filtered_vcf$genotypes, .3)

#Naive vs IFNg
covariate_names = c("sex_binary","macrophage_diff_days", "max_purity_filtered",
                    "PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))
#Test for interactions
interaction_results = testMultipleInteractions(filtered_pairs, acldl_list, filtered_vcf, formula_qtl, formula_interaction)
interaction_df = postProcessInteractionPvalues(interaction_results)
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)
write.table(interaction_hits, "results/acLDL/eQTLs/significant_interactions.txt", quote = FALSE, row.names = FALSE)

#Plot interactions
makeMultiplePlots(interaction_hits, acldl_list$cqn, filtered_vcf$genotypes, acldl_list$sample_metadata, acldl_list$gene_metadata) %>%
  savePlots("results/acLDL/eQTLs/interaction_plots/", 7,7)

#Calculate Pi1 statistic from the fastqtl results
fastqtl_min_pvalues = readRDS("results/acLDL/eQTLs/acLDL_fastqtl_min_pvalues.rds")
fastqtl_min_hits = lapply(fastqtl_min_pvalues, function(x){dplyr::filter(x, p_fdr < 0.1)})
pi1_1 = calculatePi1(fastqtl_min_pvalues$Ctrl, fastqtl_min_pvalues$AcLDL, qvalue_thresh = 0.1)
pi1_2 = calculatePi1(fastqtl_min_pvalues$AcLDL, fastqtl_min_pvalues$Ctrl, qvalue_thresh = 0.1)


