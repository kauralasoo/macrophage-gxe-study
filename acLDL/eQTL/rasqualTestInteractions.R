library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Import data
acldl_list = readRDS("results/acLDL/acLDL_combined_expression_data_covariates.rds")
acldl_list = extractConditionFromExpressionList(c("Ctrl","AcLDL"), acldl_list)

#Remove some donors
#acldl_list_filtered = acldl_list
#acldl_list_filtered$sample_metadata = acldl_list$sample_metadata %>% 
#  dplyr::filter(!(donor %in% c("voas","coio","giuo", "oefg","oarz", "hiaf","kuxp","piun", "xugn","cicb","fikt", "nusw")))
#acldl_list_filtered$cqn = acldl_list_filtered$cqn[,acldl_list_filtered$sample_metadata$sample_id]

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
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)
write.table(interaction_hits, "results/acLDL/eQTLs/significant_interactions.txt", quote = FALSE, row.names = FALSE)

#Test for interactions with filtered samples
interaction_results = testMultipleInteractions(filtered_pairs, acldl_list_filtered$cqn, acldl_list_filtered$sample_metadata, filtered_vcf, formula_qtl, formula_interaction)
interaction_df = postProcessInteractionPvalues(interaction_results)
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.1)

#Estimate parameters for top hits
interaction_params = testMultipleInteractions(dplyr::select(interaction_hits, gene_id, snp_id), acldl_list,filtered_vcf,
                         formula_qtl,formula_interaction, return = "full")
coef_matrix = lapply(interaction_params, function(x) x$interaction_model$coefficients) %>% ldply(.id = "gene_snp_pair")
write.table(coef_matrix, "results/acLDL/eQTLs/interaction_coefficients.txt", quote = FALSE, row.names = FALSE)

#Plot interactions
makeMultiplePlots(interaction_hits, acldl_list$cqn, filtered_vcf$genotypes, acldl_list$sample_metadata, acldl_list$gene_metadata) %>%
  savePlots("results/acLDL/eQTLs/interaction_plots/", 7,7)


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

