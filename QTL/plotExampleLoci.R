library("devtools")
load_all("../seqUtils/")
load_all("../wiggleplotr/")
load_all("~/software/rasqual/rasqualTools/")
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
                transcript_annotations = peak_annot, fill_palette = getGenotypePalette(), flanking_length = c(12000,300), 
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






#Analyse main SPOPL enhancer
#Plot p-values for the lead peak
peak_pvalues = tabixFetchGenesQuick("ATAC_peak_145162",tabix_file = "results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz", 
                                    gene_metadata = atac_list$gene_metadata, cis_window = 5e4)[[1]]
peak_pvals = dplyr::select(peak_pvalues, pos, p_nominal)

s_coord = 138689635 - 12000
e_coord = 138689635 + 12000
region_pvals = dplyr::filter(peak_pvals, pos > s_coord, pos < e_coord)
associated_variants = ggplot(region_pvals, aes(x = pos, y = -log(p_nominal, 10))) + 
  geom_point() + 
  theme_light() + 
  ylab("-log10(p)")
ggsave("results/ATAC/QTLs/figures/SPOPL_region_p_values.pdf", plot = associated_variants, width = 8, height = 2)

#Find peaks in the region
selected_peaks = dplyr::filter(atac_list$gene_metadata, start > s_coord, end < e_coord, chr == 2)
peak_list = list(ATAC = dplyr::transmute(selected_peaks, seqnames = chr, start, end, strand) %>% dataFrameToGRanges())
peak_annot = data_frame(transcript_id = "ATAC", gene_id = "ATAC", gene_name = "ATAC-seq", strand = "+")
flanking = c(12000, 2226)

#Make a plot of the PU.1 ChIP-Seq data
pu1_df = data_frame(sample_id = "PU1_naive", track_id = "PU.1", colour_group = "PU.1", scaling_factor = 1) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/Schmidt", paste(sample_id, ".bw", sep = "")))
spopl_region_pu1 = plotCoverage(exons = peak_list, cdss = peak_list, track_data = pu1_df, rescale_introns = FALSE, 
                                transcript_annotations = peak_annot, fill_palette = c("black"), flanking_length = flanking, 
                                connect_exons = FALSE, label_type = "peak", plot_fraction = 0.1, heights = c(0.5,0.5))
ggsave("results/ATAC/QTLs/figures/SPOPL_region_PU1.pdf", spopl_region_pu1, width = 8, height = 2)

#Construct meta-data
plotting_meta = atac_list$sample_metadata %>%
  dplyr::select(sample_id, genotype_id, condition_name) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/ATAC/", paste(sample_id, ".bw", sep = ""))) %>%
  dplyr::mutate(track_id = factor(condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))) %>%
  dplyr::left_join(sample_sizes, by = "sample_id") %>%
  dplyr::left_join(genotype_df, by = "genotype_id") %>%
  dplyr::mutate(colour_group = factor(genotype, levels = c(2,1,0)))

filtered_tracks = dplyr::filter(plotting_meta, track_id %in% c("naive"))

spopl_coverage = plotCoverage(exons = peak_list, cdss = peak_list, track_data = filtered_tracks, rescale_introns = FALSE, 
             transcript_annotations = peak_annot, fill_palette = getGenotypePalette(), flanking_length = flanking, 
             connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.5,0.5))
ggsave("results/ATAC/QTLs/figures/SPOPL_naive_coverage.pdf", spopl_coverage, width = 8, height = 3)



#Zoom out from the region a bit to include dependent peaks
selected_peaks = dplyr::filter(atac_list$gene_metadata, start > 138680000, end < 138715000, chr == 2)
peak_list = list(ATAC = dplyr::transmute(selected_peaks, seqnames = chr, start, end, strand) %>% dataFrameToGRanges())
peak_annot = data_frame(transcript_id = "ATAC", gene_id = "ATAC", gene_name = "ATAC-seq", strand = "+")
flanking = c(12000, 300)

plotting_meta = atac_list$sample_metadata %>%
  dplyr::select(sample_id, genotype_id, condition_name) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/ATAC/", paste(sample_id, ".bw", sep = ""))) %>%
  dplyr::mutate(track_id = factor(condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))) %>%
  dplyr::left_join(sample_sizes, by = "sample_id") %>%
  dplyr::left_join(genotype_df, by = "genotype_id") %>%
  dplyr::mutate(colour_group = factor(genotype, levels = c(2,1,0)))

filtered_tracks = dplyr::filter(plotting_meta, track_id %in% c("naive", "IFNg"))

spopl_coverage = plotCoverage(exons = peak_list, cdss = peak_list, track_data = filtered_tracks, rescale_introns = FALSE, 
                              transcript_annotations = peak_annot, fill_palette = getGenotypePalette(), flanking_length = flanking, 
                              connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3))
ggsave("results/ATAC/QTLs/figures/SPOPL_coverage_wider.pdf", spopl_coverage, width = 8, height = 4)


#Make new p-value plot
s_coord = 138689635 - 12000
e_coord_2 = 138712858 + 300
region_pvals = dplyr::filter(peak_pvals, pos > s_coord, pos < e_coord_2)
associated_variants = ggplot(region_pvals, aes(x = pos, y = -log(p_nominal, 10))) + 
  geom_point() + 
  theme_light() + 
  ylab("-log10(p)") + 
  scale_x_continuous(limits = c(s_coord, e_coord_2), expand = c(0,0))
ggsave("results/ATAC/QTLs/figures/SPOPL_wider_region_p_values.pdf", plot = associated_variants, width = 8, height = 2)

