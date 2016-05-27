library("devtools")
load_all("../seqUtils/")
load_all("../wiggleplotr/")
library("GenomicRanges")
library("dplyr")
library("ggplot2")

#Import data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Import genotypes
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#SPOPL enhancer
selected_peaks = dplyr::filter(atac_list$gene_metadata, gene_id %in% c("ATAC_peak_154955", "ATAC_peak_154957", "ATAC_peak_154958"))
selected_peaks = dplyr::filter(atac_list$gene_metadata, start > 138680000, end < 138715000, chr == 2)
peak_list = list(ATAC = dplyr::transmute(selected_peaks, seqnames = chr, start, end, strand) %>% dataFrameToGRanges())
peak_annot = data_frame(transcript_id = "ATAC", gene_id = "ATAC", gene_name = "ATAC-seq", strand = "+")

#Sample sizes
sample_sizes = data_frame(sample_id = colnames(atac_list$counts), scaling_factor = colSums(atac_list$counts)/1e6)

#Extract genotypes
genotype_df = vcf_file$genotypes["rs12621644",] %>% tidyVector(sample_id = "genotype_id", value_id = "genotype")

#Construct meta-data
plotting_meta = atac_list$sample_metadata %>%
  dplyr::select(sample_id, genotype_id, condition_name) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/ATAC/", paste(sample_id, ".bw", sep = ""))) %>%
  dplyr::mutate(track_id = factor(condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))) %>%
  dplyr::left_join(sample_sizes, by = "sample_id") %>%
  dplyr::left_join(genotype_df, by = "genotype_id") %>%
  dplyr::mutate(colour_group = factor(genotype, levels = c(2,1,0)))

filtered_tracks = dplyr::filter(plotting_meta, track_id %in% c("naive", "IFNg"))

#Make plot
spopl_region = plotCoverage(exons = peak_list, cdss = peak_list, track_data = filtered_tracks, rescale_introns = FALSE, 
                transcript_annotations = peak_annot, fill_palette = getGenotypePalette(), flanking_length = c(300,300), 
                connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3))
ggsave("results/ATAC/QTLs/figures/SPOPL_region.pdf", spopl_region, width = 8, height = 4)

#Make a plot of the PU.1 ChIP-Seq data
pu1_df = data_frame(sample_id = "PU1_naive", track_id = "PU.1", colour_group = "PU.1", scaling_factor = 1) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/Schmidt", paste(sample_id, ".bw", sep = "")))
spopl_region_pu1 = plotCoverage(exons = peak_list, cdss = peak_list, track_data = pu1_df, rescale_introns = FALSE, 
                            transcript_annotations = peak_annot, fill_palette = c("black"), flanking_length = c(300,300), 
                            connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.5,0.5))
ggsave("results/ATAC/QTLs/figures/SPOPL_region_PU1.pdf", spopl_region_pu1, width = 8, height = 2)







#SPOPL enhancer
selected_peaks = dplyr::filter(atac_list$gene_metadata, start > 138600000, end < 138715000, chr == 2) %>%
  dplyr::mutate(start = start - 100, end = end + 100)
peak_list = list(ATAC = dplyr::transmute(selected_peaks, seqnames = chr, start, end, strand) %>% dataFrameToGRanges())
peak_annot = data_frame(transcript_id = "ATAC", gene_id = "ATAC", gene_name = "ATAC-seq", strand = "+")

#Sample sizes
sample_sizes = data_frame(sample_id = colnames(atac_list$counts), scaling_factor = colSums(atac_list$counts)/1e6)

#Extract genotypes
genotype_df = vcf_file$genotypes["rs12621644",] %>% tidyVector(sample_id = "genotype_id", value_id = "genotype")

#Construct meta-data
plotting_meta = atac_list$sample_metadata %>%
  dplyr::select(sample_id, genotype_id, condition_name) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/ATAC/", paste(sample_id, ".bw", sep = ""))) %>%
  dplyr::mutate(track_id = factor(condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))) %>%
  dplyr::left_join(sample_sizes, by = "sample_id") %>%
  dplyr::left_join(genotype_df, by = "genotype_id") %>%
  dplyr::mutate(colour_group = factor(genotype, levels = c(2,1,0)))

filtered_tracks = dplyr::filter(plotting_meta, track_id %in% c("naive", "IFNg"))

#Make plot
spopl_region = plotCoverage(exons = peak_list, cdss = peak_list, track_data = filtered_tracks, rescale_introns = TRUE, 
                            transcript_annotations = peak_annot, fill_palette = c("#d7191c","#fdae61","#1a9641"), flanking_length = c(300,300), 
                            connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3))
ggsave("results/ATAC/QTLs/figures/SPOPL_region_wide.pdf", spopl_region, width = 8, height = 4)




