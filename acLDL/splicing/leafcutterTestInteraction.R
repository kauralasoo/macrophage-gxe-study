library("readr")
library("dplyr")
library("tidyr")
library("limma")
library("purrr")
library("devtools")
load_all("../seqUtils/")

#Import proportion data
prop_list = readRDS("results/acLDL/acLDL_combined_proportions.row_quantile.rds")
cluster_meta = dplyr::select(prop_list$gene_metadata, gene_id, cluster_id, cluster_size)

#Import min p-values
fastqtl_pvalue_list = readRDS("results/acLDL/leafcutter/leafcutter_min_pvalues.rds")

#Apply bonferroni correction for p-values within cluster
fastqtl_bonferroni = purrr::map(fastqtl_pvalue_list, ~leafcutterBonferroniCorrection(.,cluster_meta)) %>%
  purrr::map_df(identity, .id = "condition_name")

#Extract pairs
joint_pairs = dplyr::select(fastqtl_bonferroni, gene_id, snp_id) %>% unique()

#Import the VCF file
vcf_file = readRDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.rds")

#Filter VCF
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

#Set up formulas for interaction testing
covariate_names = c("norm_PC1", "norm_PC2", "norm_PC3","norm_PC4", "norm_PC5", "sex_binary")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(tbl_df(filtered_pairs), prop_list$cqn, prop_list$sample_metadata, 
                                               filtered_vcf, formula_qtl, formula_interaction)
interaction_df = postProcessInteractionPvalues(interaction_results) %>% 
  dplyr::left_join(dplyr::select(cluster_meta, gene_id, cluster_id)) %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::mutate(cluster_size = length(cluster_id)) %>% 
  dplyr::mutate(p_bonferroni = p_nominal * cluster_size) %>%
  dplyr::mutate(p_bonferroni = pmin(p_bonferroni, 1)) %>%
  dplyr::mutate(p_fdr = p.adjust(p_bonferroni, method = "fdr")) %>%
  dplyr::arrange(cluster_id, p_fdr) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(p_fdr)

#Make a couple of plots
gene_metadata = dplyr::mutate(prop_list$gene_metadata, gene_name = gene_id)
plotEQTL("17:81192593:81192736:clu_27472", "rs6565545", prop_list$tpm, vcf_file$genotypes, 
         prop_list$sample_metadata, gene_metadata)
plotEQTL("15:63504830:63532647:clu_24007", "rs62011334", prop_list$tpm, vcf_file$genotypes, 
         prop_list$sample_metadata, gene_metadata)
plotEQTL("11:110163379:110164341:clu_19512", "rs7941903", prop_list$cqn, vcf_file$genotypes, 
         prop_list$sample_metadata, gene_metadata)
plotEQTL("17:17232901:17236912:clu_26062", "rs75234140", prop_list$tpm, vcf_file$genotypes, 
         prop_list$sample_metadata, gene_metadata)
