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

#Filter transcripts for NXPH2
nxph2_tx = dplyr::filter(tx_metadata, gene_name == "NXPH2", 
                         transcript_biotype == "protein_coding", transcript_status == "KNOWN")
region_coords = c(138662000,138790000)

#Make a plot of the transcript structure
tx_plot = plotTranscripts(exons[nxph2_tx$transcript_id], cdss[nxph2_tx$transcript_id], tx_metadata, rescale_introns = FALSE, 
                          region_coords = region_coords)

#Fetch all peaks in the region and make peak annot plot
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 2, atac_list$gene_metadata)
peak_plot = plotTranscripts(peak_annot$peak_list, peak_annot$peak_list, peak_annot$peak_annot, rescale_introns = FALSE, 
                            region_coords = region_coords, connect_exons = FALSE, label_type = "peak") + dataTrackTheme()

#Make a coverage plot of the region
#Construct metadata df for wiggleplotr
atac_track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs7594476", vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive")) %>%
  dplyr::mutate(track_id = "peaks")

ATAC_coverage = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
  transcript_annotations = peak_annot$peak_annot, fill_palette = c("#525252","#525252","#525252"), 
  connect_exons = FALSE, label_type = "peak", plot_fraction = 0.1, heights = c(0.7,0.3), 
  region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "both")


#Import ATAC-seq p-values from the region
peak_pvalues = purrr::map_df(qtlResults()$atac_rasqual, ~tabixFetchGenesQuick(c("ATAC_peak_145162"),
                                 tabix_file = ., 
                                 gene_metadata = atac_list$gene_metadata, cis_window = 1e5)[[1]],
                             .id = "condition_name") %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  dplyr::arrange(p_nominal) %>%
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::mutate(track_id = factor(condition_name, levels = c("naive","IFNg")))

#Import gene p-values from the region
gene_pvalues = purrr::map_df(qtlResults()$rna_rasqual, ~tabixFetchGenesQuick(c("ENSG00000144227"),
                                          tabix_file = ., 
                                          gene_metadata = combined_expression_data$gene_metadata, cis_window = 5e5)[[1]],
                                   .id = "condition_name") %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  dplyr::arrange(p_nominal) %>%
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::mutate(track_id = factor(condition_name, levels = c("naive","IFNg")))

#Make Manhattan plots
peak_manhattan = makeManhattanPlot(peak_pvalues, region_coords, color_R2 = TRUE)
gene_manhattan = makeManhattanPlot(gene_pvalues, region_coords, color_R2 = TRUE)

#Join all plots together
joint_plot = cowplot::plot_grid(peak_manhattan, gene_manhattan, ATAC_coverage$coverage_plot, tx_plot, 
                                align = "v", ncol = 1, rel_heights = c(3,3,1.25,2.5))
ggsave("figures/main_figures/NXPH2_manhattan_plots.png", plot = joint_plot, width = 4, height = 4.5)



#### eQTL boxplot ####
#Make eQTL boxplots
NXPH2_data = constructQtlPlotDataFrame("ENSG00000144227", "rs7594476", combined_expression_data$cqn, vcf_file$genotypes, 
                                      combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  dplyr::left_join(constructGenotypeText("rs7594476", variant_information), by = "genotype_value")
NXPH2_plot = plotQtlCol(NXPH2_data)
ggsave("figures/main_figures/NXPH2_expression_boxplot.pdf", plot = NXPH2_plot, width = 2.5, height = 3.5)

#SPOPL gene
SPOPL_data = constructQtlPlotDataFrame("ENSG00000144228", "rs7594476", combined_expression_data$cqn, vcf_file$genotypes, 
                                       combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  dplyr::left_join(constructGenotypeText("rs7594476", variant_information), by = "genotype_value")
SPOPL_plot = plotQtlCol(SPOPL_data)
ggsave("figures/main_figures/SPOPL_expression_boxplot.pdf", plot = SPOPL_plot, width = 2.5, height = 3.5)



#Make a coverage plot of the whole SPOPL-NXPH2 locus
region_coords = c(138620000, 138786348)
both_tx = dplyr::filter(tx_metadata, gene_name %in% c("NXPH2","SPOPL"), 
                                    transcript_biotype == "protein_coding", transcript_status == "KNOWN")


tx_plot = plotTranscripts(exons[both_tx$transcript_id], cdss[both_tx$transcript_id], tx_metadata, rescale_introns = FALSE, region_coords = region_coords)


#Make ATAC read coverage plot of the region
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 2, atac_list$gene_metadata)

atac_track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs7594476", vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive", "IFNg"))

ATAC_coverage = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
                             transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
                             connect_exons = FALSE, label_type = "peak", plot_fraction = 0.1, heights = c(0.7,0.3), 
                             region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "both")

#Make a manhattan plot for the whole region
peak_pvals = dplyr::filter(peak_pvalues, condition_name == "naive") %>% 
  dplyr::mutate(condition_name)
peak_manhattan = makeManhattanPlot(peak_pvals, region_coords, color_R2 = TRUE)

#Make a joint plot
joint_plot = cowplot::plot_grid(peak_manhattan, ATAC_coverage$coverage_plot, tx_plot, 
                                align = "v", ncol = 1, rel_heights = c(2.2,3,2))
ggsave("figures/main_figures/NXPH2_full_region.png", joint_plot, width = 6, height = 3)



#Make a focussed coverage plot from the main peak region
region_coords = c(138688000, 138695000)
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 2, atac_list$gene_metadata)
peak_annot$peak_annot$gene_name = ""

#Make ATAC coverage plot
atac_track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs7594476", vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive", "IFNg"))
ATAC_coverage = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
                             transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
                             connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
                             region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "both")

#Make a plot of the PU.1 ChIP-Seq data
pu1_df = data_frame(sample_id = "PU1_naive", track_id = "PU.1", colour_group = "PU.1", scaling_factor = 10) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/Schmidt", paste(sample_id, ".bw", sep = "")))
spopl_region_pu1 = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = pu1_df, rescale_introns = FALSE, 
                                transcript_annotations = peak_annot$peak_annot, fill_palette = c("black"), 
                                connect_exons = FALSE, label_type = "peak", plot_fraction = 0.3, heights = c(0.5,0.5), 
                                region_coords = region_coords, return_subplots_list = TRUE)

#Make manhattan plot
pvals = makeManhattanPlot(dplyr::filter(peak_pvalues, track_id == "naive"), region_coords, color_R2 = TRUE)

joint_plot = cowplot::plot_grid(pvals, ATAC_coverage$coverage_plot, spopl_region_pu1$coverage_plot, ATAC_coverage$tx_structure, 
                                align = "v", ncol = 1, rel_heights = c(1.5,1.5,.7,1))
ggsave("figures/main_figures/NXPH2_finemapping.pdf", plot = joint_plot, width = 4, height = 3.5)





