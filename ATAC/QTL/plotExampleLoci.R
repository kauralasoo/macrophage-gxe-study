library("devtools")
load_all("../seqUtils/")
load_all("../wiggleplotr/")
load_all("~/software/rasqual/rasqualTools/")
library("GenomicRanges")
library("dplyr")
library("ggplot2")

#Import data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#SPOPL enhancer
region_coords = c(138677635, 138713158)
selected_peaks = dplyr::filter(atac_data$gene_metadata, gene_id %in% 
                                 c("ATAC_peak_154955", "ATAC_peak_154957", "ATAC_peak_154958"))
selected_peaks = dplyr::filter(atac_data$gene_metadata, start > region_coords[1], end < region_coords[2], chr == 2)
peak_annot = wiggpleplotrConstructPeakAnnotations(selected_peaks)

#Construct metadata df for wiggleplotr
meta_df = wiggleplotrConstructMetadata(atac_data$counts, atac_data$sample_metadata, "/Volumes/JetDrive/bigwigs/ATAC/")
track_data = wiggleplotrGenotypeColourGroup(meta_df, "rs12621644", vcf_file$genotypes, 1)

#Filter by condtion
filtered_tracks = dplyr::filter(track_data, track_id %in% c("naive", "IFNg"))

#Make a coverage plot of the ATAC data
spopl_region = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = filtered_tracks, rescale_introns = FALSE, 
                transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), flanking_length = c(12000,300), 
                connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
                region_coords = region_coords, return_subplots_list = TRUE)

#Make a plot of the PU.1 ChIP-Seq data
pu1_df = data_frame(sample_id = "PU1_naive", track_id = "PU.1", colour_group = "PU.1", scaling_factor = 10) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/Schmidt", paste(sample_id, ".bw", sep = "")))
spopl_region_pu1 = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = pu1_df, rescale_introns = FALSE, 
                            transcript_annotations = peak_annot$peak_annot, fill_palette = c("black"), flanking_length = c(300,300), 
                            connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.5,0.5), 
                            region_coords = region_coords, return_subplots_list = TRUE)

#Make a plot of p-values from the region
naive_peak_pvalues = tabixFetchGenesQuick(c("ATAC_peak_145162"),
                                    tabix_file = "results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz", 
                                    gene_metadata = atac_data$gene_metadata, cis_window = 1e5)
naive_peak_pvals = naive_peak_pvalues %>% 
  purrr::map_df(., ~dplyr::mutate(., track_id = "naive")) %>% 
  dplyr::select(gene_id, pos, p_nominal, track_id)

associated_variants = makeManhattanPlot(naive_peak_pvals, region_coords)
  
#Combine the two plots together
plot = cowplot::plot_grid(associated_variants, spopl_region$coverage_plot, spopl_region_pu1$coverage_plot, 
                          spopl_region$tx_structure, align = "v", rel_heights = c(0.2,0.4,0.1,0.2), ncol = 1)
ggsave("figures/main_figures/SPOPL_region_coverage.pdf", plot, width = 8, height = 8)
