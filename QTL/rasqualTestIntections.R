library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")

#Import data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")

#Import minimal p-values
min_pvalue_list = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(min_pvalue_list, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalues_df = ldply(min_pvalue_hits, .id = "condition_name")
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Calculate R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, filtered_vcf$genotypes, .8)

#Test for interactions
covariate_names = c("sex")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Run model
interaction_results = testMultipleInteractions(filtered_pairs[1:100,], atac_list, filtered_vcf, formula_qtl, formula_interaction)
interaction_df = postProcessInteractionPvalues(interaction_results)

