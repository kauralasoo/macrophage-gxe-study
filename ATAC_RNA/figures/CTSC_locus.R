library("devtools")
library("GenomicRanges")
library("dplyr")
library("ggplot2")
library("GenomicFeatures")
load_all("../seqUtils/")
load_all("../wiggleplotr/")
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
variant_information = importVariantInformation("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")


#Import transcript annotations and metadata
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_79/TranscriptDb_GRCh38_79.db")
tx_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_79/Homo_sapiens.GRCh38.79.transcript_data.rds") %>%
  dplyr::rename(transcript_id = ensembl_transcript_id,
                gene_id = ensembl_gene_id,
                gene_name = external_gene_name)
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)



##CTSC
gene_data = constructQtlPlotDataFrame("ENSG00000109861", "rs2386849", combined_expression_data$cqn, vcf_file$genotypes, 
                                      combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg_SL1344")) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::left_join(constructGenotypeText("rs2386849", variant_information), by = "genotype_value")
gene_plot = plotQTLCompact(gene_data) + ggplot2::scale_color_manual(values = conditionPalette()[c(1,4)], guide=FALSE)
ggsave("figures/main_figures/CTSC_expression_boxplot.pdf", plot = gene_plot, width = 2, height = 2.5)

peak_data = constructQtlPlotDataFrame("ATAC_peak_52208", "rs2386849", atac_list$cqn, vcf_file$genotypes, 
                                      atac_list$sample_metadata, atac_list$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg_SL1344")) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::left_join(constructGenotypeText("rs2386849", variant_information), by = "genotype_value")
peak_plot = plotQTLCompact(peak_data) + ggplot2::scale_color_manual(values = conditionPalette()[c(1,4)], guide=FALSE) + 
  ylab("Chromatin accessibility")
ggsave("figures/main_figures/CTSC_atac_boxplot.pdf", plot = peak_plot, width = 2, height = 2.5)

#Filter transcripts
#Filter transcripts for NXPH2
tx_ids = dplyr::filter(tx_metadata, gene_name == "CTSC", 
                       transcript_biotype == "protein_coding", transcript_status == "KNOWN")
region_coords = c(88293000,88370000)

tx_plot = plotTranscripts(exons["ENST00000227266"], cdss["ENST00000227266"], tx_metadata, rescale_introns = FALSE, 
                          region_coords = region_coords)

#Fetch all peaks in the region and make peak annot plot
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 11, atac_list$gene_metadata)
peak_plot = plotTranscripts(peak_annot$peak_list, peak_annot$peak_list, peak_annot$peak_annot, rescale_introns = FALSE, 
                            region_coords = region_coords, connect_exons = FALSE, label_type = "peak") + dataTrackTheme()


#Make a coverage plot of the region
#Construct metadata df for wiggleplotr
atac_track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs2386849", vcf_file$genotypes, -1) %>%
  dplyr::filter(track_id %in% c("naive","IFNg_SL1344"))%>% 
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(track_id = figure_name)

ATAC_coverage = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
                             transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(),
                             connect_exons = FALSE, label_type = "peak", plot_fraction = 0.1, heights = c(0.7,0.3), 
                             region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "both")


#Import ATAC-seq p-values from the region
peak_pvalues = purrr::map_df(qtlResults()$atac_rasqual, ~tabixFetchGenesQuick(c("ATAC_peak_52208"),
                                                                              tabix_file = ., 
                                                                              gene_metadata = atac_list$gene_metadata, cis_window = 1e5)[[1]],
                             .id = "condition_name") %>%
  dplyr::filter(condition_name %in% c("naive","IFNg_SL1344")) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(track_id = figure_name) %>%
  dplyr::arrange(p_nominal) %>%
  addR2FromLead(vcf_file$genotypes) 

#Import gene p-values from the region
gene_pvalues = purrr::map_df(qtlResults()$rna_rasqual, ~tabixFetchGenesQuick(c("ENSG00000109861"),
                                                                             tabix_file = ., 
                                                                             gene_metadata = combined_expression_data$gene_metadata, cis_window = 5e5)[[1]],
                             .id = "condition_name") %>%
  dplyr::filter(condition_name %in% c("naive","IFNg_SL1344"))  %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(track_id = figure_name) %>%
  dplyr::arrange(p_nominal) %>%
  addR2FromLead(vcf_file$genotypes)

#Make Manhattan plots
peak_manhattan = makeManhattanPlot(peak_pvalues, region_coords, color_R2 = TRUE)
gene_manhattan = makeManhattanPlot(gene_pvalues, region_coords, color_R2 = TRUE)

joint_plot = cowplot::plot_grid(peak_manhattan, gene_manhattan, ATAC_coverage$coverage_plot, tx_plot, 
                                align = "v", ncol = 1, rel_heights = c(3,3,2,2.5))
ggsave("figures/main_figures/CTSC_manhattan_plots.png", plot = joint_plot, width = 4, height = 4.5)


