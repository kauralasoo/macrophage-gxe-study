library("devtools")
library("dplyr")
library("lme4")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Import data
acldl_list = readRDS("results/acLDL/acLDL_combined_expression_data_covariates.rds")
acldl_list = extractConditionFromExpressionList(c("Ctrl","AcLDL"), acldl_list)

#Load p-values from disk
rasqual_min_pvalues = readRDS("results/acLDL/eQTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Import the VCF file
vcf_file = readRDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.rds")

#Calculate R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Filter gene-SNP pairs by R2
filtered_pairs = filterHitsR2(joint_pairs, filtered_vcf$genotypes, .8)

#Ctrl vs AcLDL
covariate_names = c("sex_binary", "PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(filtered_pairs, acldl_list$cqn, acldl_list$sample_metadata, 
                                               filtered_vcf, formula_qtl, formula_interaction, id_field_separator = "-")
interaction_df = postProcessInteractionPvalues(interaction_results, id_field_separator = "-")
saveRDS(interaction_df, "results/acLDL/eQTLs/lm_interaction_results.rds")
interaction_df = readRDS("results/acLDL/eQTLs/lm_interaction_results.rds")
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)
write.table(interaction_hits, "results/acLDL/eQTLs/significant_interactions.txt", quote = FALSE, row.names = FALSE)

#Make a Q-Q plot
qq_df = dplyr::mutate(interaction_df, p_eigen = p_nominal) %>% addExpectedPvalue()
qq_plot = ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_nominal,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")
ggsave("acLDL_figures/supplementary/Q-Q_interaction_test.pdf", plot = qq_plot, width = 5, height = 5)

#Plot interactions
makeMultiplePlots(interaction_hits, acldl_list$cqn, filtered_vcf$genotypes, acldl_list$sample_metadata, acldl_list$gene_metadata) %>%
  savePlots("results/acLDL/eQTLs/interaction_plots/", 7,7)


##### Ctrl vs AcLDL (paired design with lme4) #####
covariate_names = c("sex_binary","PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(filtered_pairs, acldl_list$cqn, acldl_list$sample_metadata, 
                                               filtered_vcf, formula_qtl, formula_interaction, id_field_separator = "-", lme4 = TRUE)
interaction_df = postProcessInteractionPvalues(interaction_results, id_field_separator = "-")
saveRDS(interaction_df, "results/acLDL/eQTLs/lme4_interaction_results.rds")
interaction_df = readRDS("results/acLDL/eQTLs/lme4_interaction_results.rds")

interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)

#Make a Q-Q plot
qq_df = dplyr::mutate(interaction_df, p_eigen = p_nominal) %>% addExpectedPvalue()
qq_plot = ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_nominal,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")
ggsave("acLDL_figures/supplementary/Q-Q_paired_interaction_test.pdf", plot = qq_plot, width = 5, height = 5)

#Plot interactions
makeMultiplePlots(interaction_hits, acldl_list$cqn, filtered_vcf$genotypes, acldl_list$sample_metadata, acldl_list$gene_metadata) %>%
  savePlots("results/acLDL/eQTLs/interaction_plots_lme4/", 7,7)

#Make fold-change plots between conditions
ctrl_samples = dplyr::filter(acldl_list$sample_metadata, condition_name == "Ctrl") %>% 
  dplyr::arrange(donor) %>% 
  dplyr::select(sample_id) %>% 
  unlist()
acldl_samples = dplyr::filter(acldl_list$sample_metadata, condition_name == "AcLDL") %>% 
  dplyr::arrange(donor) %>% 
  dplyr::select(sample_id) %>% 
  unlist()
fc_matrix = acldl_list$cqn[,acldl_samples] - acldl_list$cqn[,ctrl_samples]
sample_metadata = dplyr::filter(acldl_list$sample_metadata, condition_name == "AcLDL") %>%
  dplyr::mutate(condition_name = "AcLDL/Ctrl")

makeMultiplePlots(interaction_hits, fc_matrix, filtered_vcf$genotypes, sample_metadata, acldl_list$gene_metadata) %>%
  savePlots("results/acLDL/eQTLs/interaction_fold_change/", 7,7)

#Import fastQTL fold-change p-values from disk
FC_pvalues = importFastQTLTable("results/acLDL/fastqtl/output_FC/FC_permuted.txt.gz") %>%
  dplyr::filter(qvalue < 0.1)

makeMultiplePlots(FC_pvalues, acldl_list$cqn, vcf_file$genotypes, acldl_list$sample_metadata, acldl_list$gene_metadata) %>%
  savePlots("results/acLDL/eQTLs/FC_qtl_plots/", 7,7)



#Make an example plot of the ASE data
exon_ranges = constructExonRanges("ENSG00000141682", "rs6567134", acldl_list$gene_metadata)
sample_meta = dplyr::select(acldl_list$sample_metadata, sample_id, condition_name, genotype_id)
ase_data = fetchGeneASEData(exon_ranges, "results/acLDL/combined_ASE_counts.sorted.txt.gz", sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)


#Make plot
plotting_data = filterASEforPlotting(ase_data) %>% dplyr::filter(total_count > 10)
ggplot(plotting_data, aes(x = factor(lead_snp_value), y = ratio)) + 
  facet_wrap(~condition_name) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = .1)) +
  xlab("Feature SNP id") + 
  ylab("Reference allele ratio")


#Look for overlaps with the GWAS catalog
#Import GWAS catalog
filtered_catalog = readRDS("annotations/gwas_catalog_v1.0.1-downloaded_2016-03-02.filtered.rds")

#All GWAS overlaps
all_olaps = findGWASOverlaps(filtered_pairs, filtered_catalog, vcf_file, min_r2 = 0.6)
all_gwas_hits = dplyr::left_join(all_olaps, gene_name_map, by = "gene_id") %>%
  dplyr::select(gene_name, gene_id, snp_id, gwas_snp_id, R2, trait, gwas_pvalue)
write.table(all_gwas_hits, "results/acLDL/eQTLs/all_gwas_overlaps.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Rank traits by overlap size
ranked_traits = rankTraitsByOverlapSize(dplyr::filter(all_gwas_hits, R2 > 0.8), filtered_catalog, min_overlap = 3)

