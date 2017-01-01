library("devtools")
library("plyr")
library("dplyr")
library("ggplot2")
library("purrr")
library("GenomicFeatures")
load_all("../seqUtils/")
load_all("../wiggleplotr")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, 
                                                  levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
atac_meta_df = wiggleplotrConstructMetadata(atac_list$counts, atac_list$sample_metadata, "/Volumes/Ajamasin/bigwig/ATAC/")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import old and new variant coordinates
GRCh38_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name"))

#Import eQTL min p-values
#Load p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)

#Make a read coverage plot
region_coords = c(44175000, 44250000) 

#Construct metadata df for wiggleplotr
atac_track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs7282490", vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive"))

peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 21, atac_list$gene_metadata)

UC_coverage = plotCoverage(
  exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
  transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
  connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
  region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "both")
ggsave("figures/main_figures/UC_caQTL_overlap.png", plot = UC_coverage, width = 4, height = 4)


#Import caQTL and GWAS summary stats
#Import QTL and GWAS summary stats and convert them to the same GRCh38 coordinate space
qtl_df = data_frame(gene_id = "ATAC_peak_166661", snp_id = "rs7282490", trait = "UC")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, qtl_paths = qtlResults()$atac_fastqtl, 
                                         GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants, cis_dist = 2e5) %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::filter(condition_name %in% c("UC","naive")) %>%
  dplyr::rename(track_id = condition_name)
caqtl_manhattan_plot = makeManhattanPlot(qtl_summary, region_coords, color_R2 = TRUE, data_track = TRUE)


#Import eQTL summary stats for the ICOSLG gene
qtl_df = data_frame(gene_id = "ENSG00000160223", snp_id = "rs7282490", trait = "UC")
eqtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, qtl_paths = qtlResults()$rna_fastqtl, 
                                         GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants, cis_dist = 2e5, use_rasqual = FALSE) %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::rename(track_id = condition_name) %>%
  dplyr::mutate(track_id = "ICOSLG eQTL")

eqtl_manhattan_plot = makeManhattanPlot(eqtl_summary, region_coords, color_R2 = TRUE, data_track = TRUE)

#Make a joint plot
joint_plot = cowplot::plot_grid(caqtl_manhattan_plot, eqtl_manhattan_plot, UC_coverage$coverage_plot, UC_coverage$tx_structure, 
                                align = "v", ncol = 1, rel_heights = c(4,2,2,1))
ggsave("figures/main_figures/UC_caQTL_overlap.png", plot = joint_plot, width = 4, height = 4)



#Make a coverage plot from the narrow region around the peak
region_coords = c(44196000, 44198000)

#Construct metadata df for wiggleplotr
atac_track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs7282490", vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive"))

peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 21, atac_list$gene_metadata)

UC_coverage1 = plotCoverage(
  exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
  transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
  connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
  region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "both")

caqtl_manhattan_plot = makeManhattanPlot(qtl_summary, region_coords, color_R2 = TRUE, data_track = TRUE)
narrow_joint_plot = cowplot::plot_grid(caqtl_manhattan_plot, UC_coverage1$coverage_plot, UC_coverage1$tx_structure, 
                   align = "v", ncol = 1, rel_heights = c(3,3,1))
ggsave("figures/supplementary/UC_caQTL_overlap_narrow.pdf", plot = narrow_joint_plot, width = 4, height = 4)


#Check if the causal SNPs disrupt any known TF motifs
#Import motif matches
motif_metadata = readRDS("results/ATAC/cisBP/cisBP_motif_metadata.rds") %>%
  dplyr::transmute(motif_id = Motif_ID, tf_name = TF_Name, tf_count = TF_count)
motif_disruptions = importMotifDisruptions("results/ATAC/motif_analysis/motif_disruption.txt") %>%
  dplyr::left_join(motif_metadata, by = "motif_id")

#Filter by SNP ID
motif_hits = dplyr::filter(motif_disruptions, snp_id %in% causal_variants) %>%
  dplyr::filter(max_rel_score > 0.8) %>% dplyr::arrange(-abs(rel_diff))




#Make a QTL plot for ICOSLG
#Make a boxplot of the eQTL data
ICOSLG_data = constructQtlPlotDataFrame("ENSG00000160223", "rs4819387", 
                                       combined_expression_data$cqn, 
                                       vcf_file$genotypes, 
                                       combined_expression_data$sample_metadata, 
                                       combined_expression_data$gene_metadata) %>%
  dplyr::left_join(constructGenotypeText("rs4819387", GRCh38_variants), by = "genotype_value")
plotQtlCol(ICOSLG_data, scales = "free_y")

#Calculate R2 between the GWAS lead variant and both eQTL lead variant
calculatePairR2("rs4819387", "rs7282490", vcf_file$genotypes)
calculatePairR2("rs6518352", "rs7282490", vcf_file$genotypes)





#Import caQTL and GWAS summary stats
#Import QTL and GWAS summary stats and convert them to the same GRCh38 coordinate space
qtl_df = data_frame(gene_id = "ENSG00000160223", snp_id = "rs7282490", trait = "UC")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, qtl_paths = qtlResults()$rna_fastqtl, 
                                         GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants, cis_dist = 2e5, use_rasqual = FALSE) %>%
  arrange(p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::rename(track_id = condition_name)
region_coords = c(min(qtl_summary$pos), max(qtl_summary$pos))
region_coords = c(44175000, 44255000) 
makeManhattanPlot(qtl_summary, region_coords, color_R2 = TRUE, data_track = TRUE)



#Mark the likely causal variants
causal_variants = c("rs4819387","rs6518352", "rs7282490")
names(causal_variants) = c("c1","c2")

#Extract relevant genotypes from vcf
extracted_genotypes = t(vcf_file$genotypes[causal_variants,]) %>% data.frame()
extracted_genotypes = dplyr::mutate(extracted_genotypes, genotype_id = rownames(extracted_genotypes)) %>%
  dplyr::select(genotype_id, everything())


#Make initial formula
covariate_names = c("PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6", "sex_binary")
formula = as.formula(paste("norm_exp ~ rs4819387 ",
                           paste(covariate_names, collapse = " + "), sep = "+ "))
formula2 = as.formula(paste("norm_exp ~ rs4819387 + rs6518352 ",
                           paste(covariate_names, collapse = " + "), sep = "+ "))
formula3 = as.formula(paste("norm_exp ~ rs4819387 + rs7282490 ",
                            paste(covariate_names, collapse = " + "), sep = "+ "))

#Make dataset
data = dplyr::left_join(ICOSLG_data, extracted_genotypes, by = "genotype_id") %>%
  dplyr::filter(condition_name == "SL1344")

#Get fits
fit1 = lm(formula, data)
fit2 = lm(formula2, data)
fit3 = lm(formula3, data)

data = dplyr::mutate(data, resid = as.numeric(lmfit$residuals))


summary(lm(resid ~ rs28834970, data))






#Import minimal p-values
rasqual_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)


#Load p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::left_join(gene_name_map)

dplyr::filter(min_pvalues_df, gene_name == "ICOSLG")
cor(vcf_file$genotypes["rs4819387",], vcf_file$genotypes["rs7282490",])^2
cor(vcf_file$genotypes["rs4511811",], vcf_file$genotypes["rs7282490",])^2

dplyr::filter(min_pvalues_df, gene_name == "TRAPPC10")


