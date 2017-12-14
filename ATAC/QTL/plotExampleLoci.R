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
atac_meta_df = wiggleplotrConstructMetadata(atac_data$counts, atac_data$sample_metadata, "/Volumes/JetDrive/bigwigs/ATAC/")
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#SPOPL enhancer
region_coords = c(138677635, 138713158)
selected_peaks = dplyr::filter(atac_data$gene_metadata, gene_id %in% 
                                 c("ATAC_peak_154955", "ATAC_peak_154957", "ATAC_peak_154958"))
selected_peaks = dplyr::filter(atac_data$gene_metadata, start > region_coords[1], end < region_coords[2], chr == 2)
peak_annot = wiggpleplotrConstructPeakAnnotations(selected_peaks)

#Check if the causal SNPs disrupt any known TF motifs
#Import motif matches
motif_metadata = readRDS("results/ATAC/cisBP/cisBP_motif_metadata.rds") %>%
  dplyr::transmute(motif_id = Motif_ID, tf_name = TF_Name, tf_count = TF_count)
motif_disruptions = importMotifDisruptions("results/ATAC/motif_analysis/motif_disruption.txt") %>%
  dplyr::left_join(motif_metadata, by = "motif_id")

#Filter by SNP ID
motif_hits = dplyr::filter(motif_disruptions, snp_id %in% c("rs7594476")) %>%
  dplyr::filter(max_rel_score > 0.8) %>% dplyr::arrange(-abs(rel_diff))


#Construct metadata df for wiggleplotr
track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs7594476", vcf_file$genotypes, 1)

#Filter by condtion
filtered_tracks = dplyr::filter(track_data, track_id %in% c("naive", "IFNg"))

#Make a coverage plot of the ATAC data
spopl_region = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = filtered_tracks, rescale_introns = FALSE, 
                transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
                connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
                region_coords = region_coords, return_subplots_list = TRUE)

#Make a plot of the PU.1 ChIP-Seq data
pu1_df = data_frame(sample_id = "PU1_naive", track_id = "PU.1", colour_group = "PU.1", scaling_factor = 10) %>%
  dplyr::mutate(bigWig = file.path("/Volumes/JetDrive/bigwigs/Schmidt", paste(sample_id, ".bw", sep = "")))
spopl_region_pu1 = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = pu1_df, rescale_introns = FALSE, 
                            transcript_annotations = peak_annot$peak_annot, fill_palette = c("black"), 
                            connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.5,0.5), 
                            region_coords = region_coords, return_subplots_list = TRUE)

#Make a plot of p-values from the region
naive_peak_pvalues = tabixFetchGenesQuick(c("ATAC_peak_145162"),
                                    tabix_file = qtlResults()$atac_rasqual$naive, 
                                    gene_metadata = atac_data$gene_metadata, cis_window = 1e5)
naive_peak_pvals = naive_peak_pvalues %>% 
  purrr::map_df(., ~dplyr::mutate(., track_id = "naive")) %>% 
  dplyr::select(gene_id, pos, p_nominal, track_id)

associated_variants = makeManhattanPlot(naive_peak_pvals, region_coords)
  
#Combine the two plots together
plot = cowplot::plot_grid(associated_variants, spopl_region$coverage_plot, spopl_region_pu1$coverage_plot, 
                          spopl_region$tx_structure, align = "v", rel_heights = c(0.2,0.4,0.1,0.2), ncol = 1)
ggsave("figures/main_figures/SPOPL_region_coverage.pdf", plot, width = 4.5, height = 5)



#Make Manhattan plots of the SPOPL region
#Import chromatin QTL pvalues
region_coords = c(138643635, 138743158)
naive_peak_pvalues = tabixFetchGenesQuick(c("ATAC_peak_145162"),
                                          tabix_file = qtlResults()$atac_rasqual$naive, 
                                          gene_metadata = atac_data$gene_metadata, cis_window = 1e5) %>%
  purrr::map_df(., ~dplyr::mutate(., track_id = "naive caQTL"))
IFNg_peak_pvalues = tabixFetchGenesQuick(c("ATAC_peak_145162"),
                                          tabix_file = qtlResults()$atac_rasqual$IFNg, 
                                          gene_metadata = atac_data$gene_metadata, cis_window = 1e5) %>%
  purrr::map_df(., ~dplyr::mutate(., track_id = "IFNg caQTL"))

#Import gene expression p-values
naive_gene_pvalues = tabixFetchGenesQuick(c("ENSG00000144228"),
                                          tabix_file = qtlResults()$rna_rasqual$naive, 
                                          gene_metadata = combined_expression_data$gene_metadata, cis_window = 5e5) %>%
  purrr::map_df(., ~dplyr::mutate(., track_id = "naive eQTL"))
IFNg_gene_pvalues = tabixFetchGenesQuick(c("ENSG00000144228"),
                                         tabix_file = qtlResults()$rna_rasqual$IFNg, 
                                         gene_metadata = combined_expression_data$gene_metadata, cis_window = 5e5) %>%
  purrr::map_df(., ~dplyr::mutate(., track_id = "IFNg eQTL"))

peak_pvals = dplyr::bind_rows(naive_peak_pvalues, IFNg_peak_pvalues, naive_gene_pvalues, IFNg_gene_pvalues) %>% 
  dplyr::select(gene_id, pos, p_nominal, track_id) %>%
  dplyr::mutate(track_id = factor(track_id, levels = c("naive caQTL", "naive eQTL","IFNg caQTL","IFNg eQTL")))

spopl_manhattan = makeManhattanPlot(peak_pvals, region_coords)
ggsave("figures/main_figures/SPOPL_manhattan.pdf", spopl_manhattan, width = 4, height = 6)



#Some simple QTL plots
plot = plotEQTL("ENSG00000150782", "rs71478720", combined_expression_data$cqn, vcf_file$genotypes, 
                combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
ggsave("IL18.pdf", plot = plot, width = 8, height = 8)


plot = plotEQTL("ENSG00000150782", "rs10891343", combined_expression_data$cqn, vcf_file$genotypes, 
                combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
ggsave("IL18_lead.pdf", plot = plot, width = 8, height = 8)


plot = plotEQTL("ENSG00000091106", "rs385076", combined_expression_data$cqn, vcf_file$genotypes, 
                combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
ggsave("NLRC4.pdf", plot = plot, width = 8, height = 8)

plot = plotEQTL("ENSG00000150782", "rs385076", combined_expression_data$cqn, vcf_file$genotypes, 
                combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
ggsave("IL18_rs385076.pdf", plot = plot, width = 8, height = 8)

