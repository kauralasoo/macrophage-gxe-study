library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
library("ggplot2")
library("wiggleplotr")
library("GenomicFeatures")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("../wiggleplotr/")


#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, 
                                                  levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
atac_meta_df = wiggleplotrConstructMetadata(atac_list$counts, atac_list$sample_metadata, "/Volumes/Ajamasin/bigwig/ATAC/")

#Import old and new variant coordinates
GRCh38_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name"))

#Import transcript annotations and metadata
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_79/TranscriptDb_GRCh38_79.db")
tx_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_79/Homo_sapiens.GRCh38.79.transcript_data.rds") %>%
  dplyr::rename(transcript_id = ensembl_transcript_id,
                gene_id = ensembl_gene_id,
                gene_name = external_gene_name)
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

#Import rescued pair
filtered_overlaps = read.table("results/ATAC_RNA_overlaps/caQTL_eQTL_rescued_pairs.txt", stringsAsFactors = FALSE, header = TRUE)

#Export CTSB
#Import caQTL and GWAS summary stats
#Import QTL and GWAS summary stats and convert them to the same GRCh38 coordinate space
qtl_df = data_frame(phenotype_id = "ATAC_peak_260182", snp_id = "rs11997338", trait = "SLE")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, qtl_paths = qtlResults()$atac_fastqtl, 
                                         GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants, 
                                         cis_dist = 2e5, gwas_dir = "~/datasets/Inflammatory_GWAS/") %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::filter(condition_name %in% c("SLE","naive","IFNg_SL1344")) %>%
  dplyr::mutate(track_id = as.character(condition_name)) %>%
  dplyr::mutate(track_id = ifelse(track_id == "naive", "N", ifelse(track_id == "IFNg_SL1344", "I+S",track_id))) %>%
  dplyr::mutate(track_id = factor(track_id, levels = c("SLE","N", "I+S")))
region_coords = c(min(qtl_summary$pos), max(qtl_summary$pos))
caqtl_manhattan_plot = makeManhattanPlot(qtl_summary, region_coords, color_R2 = TRUE, data_track = TRUE)
ggsave("figures/supplementary/CTSB_caQTL_fastqtl.pdf", caqtl_manhattan_plot, width = 6, height = 8)


#Import eQTL summaries (FastQTL)
qtl_df = data_frame(phenotype_id = "ENSG00000164733", snp_id = "rs11997338", trait = "SLE")
eqtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, qtl_paths = qtlResults()$rna_fastqtl, 
                                          GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants,
                                          cis_dist = 2e5, gwas_dir = "~/datasets/Inflammatory_GWAS/") %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::filter(condition_name %in% c("SLE","naive","IFNg_SL1344")) %>%
  dplyr::mutate(track_id = as.character(condition_name)) %>%
  dplyr::mutate(track_id = ifelse(track_id == "naive", "N", ifelse(track_id == "IFNg_SL1344", "I+S",track_id))) %>%
  dplyr::mutate(track_id = factor(track_id, levels = c("SLE","N", "I+S")))
region_coords = c(min(qtl_summary$pos), max(qtl_summary$pos))
eqtl_manhattan_plot = makeManhattanPlot(eqtl_summary, region_coords, color_R2 = TRUE, data_track = TRUE)
ggsave("figures/supplementary/CTSB_eQTL_fastqtl.pdf", eqtl_manhattan_plot, width = 6, height = 8)


#Import eQTL summaries (RASQUAL)
qtl_df = data_frame(phenotype_id = "ENSG00000164733", snp_id = "rs11997338", trait = "SLE")
eqtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, qtl_paths = qtlResults()$rna_rasqual, 
                                          GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants,
                                          cis_dist = 2e5, gwas_dir = "~/datasets/Inflammatory_GWAS/", type = "RASQUAL") %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::mutate(track_id = condition_name)
region_coords = c(min(qtl_summary$pos), max(qtl_summary$pos))
gwas_manhattan_plot = makeManhattanPlot(dplyr::filter(eqtl_summary, condition_name == "SLE"), 
                                        region_coords, color_R2 = TRUE, data_track = TRUE)
eqtl_manhattan_plot = makeManhattanPlot(dplyr::filter(eqtl_summary, condition_name != "SLE"), 
                                        region_coords, color_R2 = TRUE,data_track = TRUE)

#Join all plots together
joint_plot = cowplot::plot_grid(gwas_manhattan_plot, eqtl_manhattan_plot, 
                                align = "v", ncol = 1, rel_heights = c(1,4))
ggsave("figures/supplementary/CTSB_eQTL_rasqual.pdf", joint_plot, width = 6, height = 8)



#Make QTL boxplots
gene_data = constructQtlPlotDataFrame("ENSG00000164733", "rs11997338", combined_expression_data$cqn, vcf_file$genotypes, 
                                      combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg_SL1344")) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::left_join(constructGenotypeText("rs11997338", GRCh38_variants), by = "genotype_value")
gene_plot = plotQTLCompact(gene_data) + ggplot2::scale_color_manual(values = conditionPalette()[c(1,4)], guide=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(paste(italic("CTSB"), " expression")))
ggsave("figures/main_figures//CTSB_expression_boxplot.pdf", plot = gene_plot, width = 2, height = 2.5)

peak_data = constructQtlPlotDataFrame("ATAC_peak_260182", "rs11997338", atac_list$cqn, vcf_file$genotypes, 
                                      atac_list$sample_metadata, atac_list$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg_SL1344")) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::left_join(constructGenotypeText("rs11997338", GRCh38_variants), by = "genotype_value")
peak_plot = plotQTLCompact(peak_data) + ggplot2::scale_color_manual(values = conditionPalette()[c(1,4)], guide=FALSE) + 
  ylab("Chromatin accessibility") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("figures/main_figures/CTSB_chromatin_boxplot.pdf", plot = peak_plot, width = 2, height = 2.5)


#Redo coloc analysis
eqtl_df = data_frame(phenotype_id = "ENSG00000164733", snp_id = "rs11997338")
eqtl_coloc = colocMolecularQTLs(eqtl_df,qtl_summary_path = qtlResults()$rna_fastqtl$IFNg_SL1344, 
              gwas_summary_path = "../../datasets/Inflammatory_GWAS/Systemic_lupus_erythematosus_Bentham_2015_NatGen_GWAS.sorted.txt.gz",
              GRCh37_variants, GRCh38_variants, QTLTools = FALSE, cis_dist = 2e5)


eqtl_df = data_frame(phenotype_id = "ENSG00000164733", snp_id = "rs7835672")
eqtl_coloc = colocMolecularQTLs(eqtl_df,qtl_summary_path = qtlResults()$rna_fastqtl$IFNg_SL1344, 
                                gwas_summary_path = "../../datasets/Inflammatory_GWAS/Systemic_lupus_erythematosus_Bentham_2015_NatGen_GWAS.sorted.txt.gz",
                                GRCh37_variants, GRCh38_variants, QTLTools = FALSE, cis_dist = 2e5)


#Filter transcripts
#Filter transcripts for NXPH2
tx_ids = dplyr::filter(tx_metadata, gene_name == "CTSB", 
                       transcript_biotype == "protein_coding", transcript_status == "KNOWN")
region_coords = c(11842520,11893000)

tx_plot = plotTranscripts(exons["ENST00000453527"], cdss["ENST00000453527"], tx_metadata, rescale_introns = FALSE, 
                          region_coords = region_coords)

#Fetch all peaks in the region and make peak annot plot
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 8, atac_list$gene_metadata)
peak_plot = plotTranscripts(peak_annot$peak_list, peak_annot$peak_list, peak_annot$peak_annot, rescale_introns = FALSE, 
                            region_coords = region_coords, connect_exons = FALSE, transcript_label = FALSE) + dataTrackTheme()


#Make a coverage plot of the region
#Construct metadata df for wiggleplotr
atac_track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs11997338", vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive","IFNg_SL1344"))%>% 
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(track_id = figure_name)

ATAC_coverage = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
                             transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(),
                             connect_exons = FALSE, transcript_label = FALSE, plot_fraction = 0.1, heights = c(0.7,0.3), 
                             region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "both")
saveRDS(ATAC_coverage, "figures/main_figures/CTSB_coverage_data.rds")


#Make manhattan plots
caQTL_plot = makeManhattanPlot(qtl_summary, region_coords, color_R2 = TRUE, data_track = TRUE)
eqtl_plot = makeManhattanPlot(dplyr::filter(eqtl_summary, track_id != "SLE"), region_coords, color_R2 = TRUE, data_track = TRUE)

#Make a joint plot
joint_plot = cowplot::plot_grid(caQTL_plot, eqtl_plot, ATAC_coverage$coverage_plot, tx_plot, 
                                align = "v", ncol = 1, rel_heights = c(3,2,2,1.5))
ggsave("figures/main_figures/CTSB_joint_plot.pdf", plot = joint_plot, width = 4, height = 7)
