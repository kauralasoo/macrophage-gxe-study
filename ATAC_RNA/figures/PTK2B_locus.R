library("devtools")
library("GenomicRanges")
library("dplyr")
library("ggplot2")
library("GenomicFeatures")
load_all("../seqUtils/")
load_all("../wiggleplotr/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
atac_meta_df = wiggleplotrConstructMetadata(atac_data$counts, atac_data$sample_metadata, "/Volumes/Ajamasin/bigwig/ATAC/")
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")
variant_information = importVariantInformation("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")


#### PTK2B region ####
ptk2b_coords = c(27331000, 27345000)
selected_peaks = dplyr::filter(atac_data$gene_metadata, start > ptk2b_coords[1], end < ptk2b_coords[2], chr == 8)
peak_annot = wiggpleplotrConstructPeakAnnotations(selected_peaks)

#Add lead SNP genotype
track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs28834970", vcf_file$genotypes, -1) %>%
  dplyr::filter(track_id %in% c("naive")) %>%
  dplyr::mutate(track_id = "N")

#Make fragment coverage plot
ptk2b_plot = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = track_data, rescale_introns = FALSE, 
                          transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
                          connect_exons = FALSE, transcript_label = FALSE, plot_fraction = 0.2, heights = c(0.7,0.3), 
                          region_coords = ptk2b_coords, return_subplots_list = TRUE)
atac_coverage = ptk2b_plot$coverage_plot + coord_cartesian(ylim = c(0,4))

#Make a plot of the CEBPb ChIP-Seq data
cebpb_df = data_frame(sample_id = "CEBPbeta_ctrl_204", track_id = "CEBPb", colour_group = "CEBPb", scaling_factor = 10) %>%
  dplyr::mutate(bigWig = file.path("~/datasets/bigwigs/", paste(sample_id, ".bw", sep = "")))

cebpb_plot = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = cebpb_df, rescale_introns = FALSE, 
                          transcript_annotations = peak_annot$peak_annot, fill_palette = c("black"), 
                          connect_exons = FALSE,  transcript_label = FALSE, plot_fraction = 0.2, heights = c(0.7,0.3), 
                          region_coords = ptk2b_coords, return_subplots_list = TRUE)
cebpb_coverage = cebpb_plot$coverage_plot + coord_cartesian(ylim = c(0,20))

#Make a plot of p-values from the region
naive_peak_pvalues = tabixFetchGenesQuick(c("ATAC_peak_261927"),
                                          tabix_file = qtlResults()$atac_rasqual$naive, 
                                          gene_metadata = atac_data$gene_metadata, cis_window = 1e5)
naive_peak_pvals = naive_peak_pvalues %>% 
  purrr::map_df(., ~dplyr::mutate(., track_id = "N")) %>% 
  dplyr::select(gene_id, pos, p_nominal, track_id)
associated_variants = makeManhattanPlot(naive_peak_pvals, ptk2b_coords)

#Join all of the subplots together
ptk2b_full = cowplot::plot_grid(associated_variants, atac_coverage, cebpb_coverage, 
                                ptk2b_plot$tx_structure, align = "v", rel_heights = c(0.25,0.25,0.1,0.2), ncol = 1)
ggsave("figures/supplementary/PTK2B_region_coverage.pdf", ptk2b_full, width = 4.5, height = 5)




#Make a boxplot of PTK2B expression
gene_data = constructQtlPlotDataFrame("ENSG00000120899", "rs28834970", combined_expression_data$cqn, vcf_file$genotypes, 
                                      combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::left_join(constructGenotypeText("rs28834970", variant_information), by = "genotype_value")
gene_plot = plotQTLCompact(gene_data) + 
  ggplot2::scale_color_manual(values = conditionPalette()[c(1,4)], guide=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(paste(italic("PTK2B"), " expression")))
ggsave("figures/supplementary/PTK2B_boxplot.pdf", gene_plot, width = 3, height = 3)


#Make a read coverage plot for RNA-Seq
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_79/TranscriptDb_GRCh38_79.db")
tx_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_79/Homo_sapiens.GRCh38.79.transcript_data.rds") %>%
  dplyr::rename(transcript_id = ensembl_transcript_id,
                gene_id = ensembl_gene_id,
                gene_name = external_gene_name)
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

#Filter transcripts for PTK2B
ptk2b_tx = dplyr::filter(tx_metadata, gene_name == "PTK2B", 
                         transcript_biotype == "protein_coding", transcript_status == "KNOWN")
plotTranscripts(exons[ptk2b_tx$transcript_id], cdss[ptk2b_tx$transcript_id], tx_metadata, rescale_introns = FALSE)



