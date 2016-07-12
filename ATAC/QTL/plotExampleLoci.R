library("devtools")
load_all("../seqUtils/")
load_all("../wiggleplotr/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")
library("GenomicRanges")
library("dplyr")
library("ggplot2")


#Import data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
atac_meta_df = wiggleplotrConstructMetadata(atac_data$counts, atac_data$sample_metadata, "/Volumes/JetDrive/bigwigs/ATAC/")

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")


#SPOPL enhancer
region_coords = c(138677635, 138713158)
selected_peaks = dplyr::filter(atac_data$gene_metadata, gene_id %in% 
                                 c("ATAC_peak_154955", "ATAC_peak_154957", "ATAC_peak_154958"))
selected_peaks = dplyr::filter(atac_data$gene_metadata, start > region_coords[1], end < region_coords[2], chr == 2)
peak_annot = wiggpleplotrConstructPeakAnnotations(selected_peaks)

#Construct metadata df for wiggleplotr
track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs12621644", vcf_file$genotypes, 1)

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





#### PTK2B region ####
ptk2b_coords = c(27331000, 27345000)
selected_peaks = dplyr::filter(atac_data$gene_metadata, start > ptk2b_coords[1], end < ptk2b_coords[2], chr == 8)
peak_annot = wiggpleplotrConstructPeakAnnotations(selected_peaks)

#Add lead SNP genotype
track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs28834970", vcf_file$genotypes, -1) %>%
  dplyr::filter(track_id %in% c("naive", "IFNg"))

#Make fragment coverage plot
ptk2b_plot = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = track_data, rescale_introns = FALSE, 
                            transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
                            connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
                            region_coords = ptk2b_coords, return_subplots_list = TRUE)

#Make a plot of the CEBPb ChIP-Seq data
cebpb_df = data_frame(sample_id = "CEBPbeta_ctrl_204", track_id = "CEBPb", colour_group = "CEBPb", scaling_factor = 10) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/Reschen/", paste(sample_id, ".bw", sep = "")))

cebpb_plot = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = cebpb_df, rescale_introns = FALSE, 
                            transcript_annotations = peak_annot$peak_annot, fill_palette = c("black"), 
                            connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
                            region_coords = ptk2b_coords, return_subplots_list = TRUE)

#Make a plot of p-values from the region
naive_peak_pvalues = tabixFetchGenesQuick(c("ATAC_peak_261927"),
                                          tabix_file = qtlResults()$atac_rasqual$naive, 
                                          gene_metadata = atac_data$gene_metadata, cis_window = 1e5)
naive_peak_pvals = naive_peak_pvalues %>% 
  purrr::map_df(., ~dplyr::mutate(., track_id = "naive")) %>% 
  dplyr::select(gene_id, pos, p_nominal, track_id)
associated_variants = makeManhattanPlot(naive_peak_pvals, ptk2b_coords)

#Join all of the subplots together
ptk2b_full = cowplot::plot_grid(associated_variants, ptk2b_plot$coverage_plot, cebpb_plot$coverage_plot, 
                                ptk2b_plot$tx_structure, align = "v", rel_heights = c(0.1,0.4,0.1,0.2), ncol = 1)
ggsave("figures/main_figures/PTK2B_region_coverage.pdf", ptk2b_full, width = 8, height = 8)


#### PTK2B second QTL region ####
ptk2b_coords = c(27378000, 27384000)
selected_peaks = dplyr::filter(atac_data$gene_metadata, start > ptk2b_coords[1], end < ptk2b_coords[2], chr == 8)
peak_annot = wiggpleplotrConstructPeakAnnotations(selected_peaks)

#Add lead SNP genotype
track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs10086852", vcf_file$genotypes, 1)

#Make fragment coverage plot
ptk2b_plot = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = track_data, rescale_introns = FALSE, 
                          transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
                          connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
                          region_coords = ptk2b_coords, return_subplots_list = TRUE)

#Make a plot of p-values from the region
IFNg_SL1344_peak_pvalues = tabixFetchGenesQuick(c("ATAC_peak_261942"),
                                          tabix_file = qtlResults()$atac_rasqual$IFNg_SL1344, 
                                          gene_metadata = atac_data$gene_metadata, cis_window = 1e5)
IFNg_SL1344_peak_pvals = IFNg_SL1344_peak_pvalues %>% 
  purrr::map_df(., ~dplyr::mutate(., track_id = "naive")) %>% 
  dplyr::select(gene_id, pos, p_nominal, track_id)
enh2_associated_variants = makeManhattanPlot(IFNg_SL1344_peak_pvals, ptk2b_coords)

#Combine plots together
ptk2b_enh2 = cowplot::plot_grid(enh2_associated_variants, ptk2b_plot$coverage_plot, 
                   ptk2b_plot$tx_structure, align = "v", rel_heights = c(0.1,0.4,0.2), ncol = 1)
ggsave("figures/main_figures/PTK2B_second_enhancer.pdf", ptk2b_enh2, width = 8, height = 8)

