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
region_coords = c(138670000,138782000)

#Make a plot of the transcript structure
tx_plot = plotTranscripts(exons[nxph2_tx$transcript_id], cdss[nxph2_tx$transcript_id], tx_metadata, rescale_introns = FALSE, 
                          region_coords = region_coords)

#Fetch all peaks in the region and make peak annot plot
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 2, atac_list$gene_metadata)
peak_plot = plotTranscripts(peak_annot$peak_list, peak_annot$peak_list, peak_annot$peak_annot, rescale_introns = FALSE, 
                            region_coords = region_coords, connect_exons = FALSE, label_type = "peak") + dataTrackTheme()


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
joint_plot = cowplot::plot_grid(peak_manhattan, gene_manhattan, peak_plot, tx_plot, 
                                align = "v", ncol = 1, rel_heights = c(3,3,1.25,2.5))
ggsave("figures/main_figures/NXPH2_manhattan_plots.pdf", plot = joint_plot, width = 4, height = 4.5)



#### eQTL boxplot ####
#Make eQTL boxplot
NXPH2_data = constructQtlPlotDataFrame("ENSG00000144227", "rs7594476", combined_expression_data$cqn, vcf_file$genotypes, 
                                      combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg"))
NXPH2_plot = plotQtlCol(NXPH2_data)
ggsave("figures/main_figures/NXPH2_expression_boxplot.pdf", plot = NXPH2_plot, width = 2.5, height = 3.5)


