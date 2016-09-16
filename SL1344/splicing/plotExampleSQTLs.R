library("GenomicFeatures")
library("dplyr")
library("devtools")
load_all("../reviseAnnotations/")
load_all("../wiggleplotr/")
load_all("../seqUtils/")
load_all("../../software/rasqual/rasqualTools/")
library("ggplot2")
library("SummarizedExperiment")

#Import transcript annotations
gene_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_85/Homo_sapiens.GRCh38.85.compiled_tx_metadata.filtered.rds")  
txdb = AnnotationDbi::loadDb("../../annotations/GRCh38/genes/Ensembl_85/TranscriptDb_GRCh38_85.db")
exons = GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)

#Filter gene annotations
filtered_metadata = dplyr::filter(gene_metadata, transcript_biotype %in% c("lincRNA", "protein_coding")) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(n_transcripts = length(ensembl_gene_id)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(n_transcripts > 1) %>%
  reviseAnnotations::markLongestTranscripts()

#Make some plots
plotting_annotations = dplyr::select(filtered_metadata, ensembl_transcript_id, ensembl_gene_id, external_gene_name, strand) %>% 
  dplyr::rename(transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import gene level data
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Make coverage metadata
#Import data
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
rna_meta_str1_df = wiggleplotrConstructMetadata(combined_expression_data$counts, 
            combined_expression_data$sample_metadata, "/Volumes/Ajamasin/bigwig/RNA/", bigWig_suffix = ".str1.bw") %>%
  dplyr::mutate(scaling_factor = 1)
rna_meta_str2_df = wiggleplotrConstructMetadata(combined_expression_data$counts, 
            combined_expression_data$sample_metadata, "/Volumes/Ajamasin/bigwig/RNA/", bigWig_suffix = ".str2.bw")

##### Plot IRF5 sQTLs ####

#Extract gene data from annotations
gene_data = reviseAnnotations::extractGeneData("ENSG00000128604", filtered_metadata, exons, cdss)

#Extend truncated transcripts until the longest transcript
gene_extended_tx = reviseAnnotations::extendTranscriptsPerGene(gene_data$metadata, gene_data$exons, gene_data$cdss)
gene_data_ext = reviseAnnotations::replaceExtendedTranscripts(gene_data, gene_extended_tx)

#Plot raw protein coding transcript annotations
irf5_raw_transcripts = wiggleplotr::plotTranscripts(gene_data$exons, gene_data$cdss, 
                                                    plotting_annotations, rescale_introns = TRUE)
ggsave("figures/supplementary/IRF5_raw_transcripts.pdf", plot = irf5_raw_transcripts, width = 7, height = 5)

irf5_extended_transcripts = wiggleplotr::plotTranscripts(gene_data_ext$exons, 
                          gene_data_ext$cdss, plotting_annotations, rescale_introns = TRUE)
ggsave("figures/supplementary/IRF5_extended_transcripts.pdf", plot = irf5_extended_transcripts, width = 7, height = 5)

#Plot events constructed from the IRF5 transcripts
alt_events = constructAlternativeEvents(gene_data_ext$exons, "ENSG00000128604")
up_events = wiggleplotr::plotTranscripts(alt_events[[1]]$upstream, alt_events[[1]]$upstream, plotting_annotations, rescale_introns = TRUE)
down_events = wiggleplotr::plotTranscripts(alt_events[[1]]$downstream, alt_events[[1]]$downstream, plotting_annotations, rescale_introns = TRUE)
contained_events = wiggleplotr::plotTranscripts(alt_events[[1]]$contained, alt_events[[1]]$contained, plotting_annotations, rescale_introns = TRUE)
ggsave("figures/supplementary/IRF5_upstream.pdf", plot = up_events, width = 7, height = 4)
ggsave("figures/supplementary/IRF5_downstream.pdf", plot = down_events, width = 7, height = 3)
ggsave("figures/supplementary/IRF5_contained.pdf", plot = contained_events, width = 7, height = 2)



#### Make coverage plots ####
#Calculate union exons
union_exons = list("Union" = listUnion(gene_data_ext$exons))
union_cds = list("Union" = listUnion(gene_data_ext$cdss))
union_annotations = dplyr::bind_rows(plotting_annotations, 
  data_frame(transcript_id = "Union", gene_id = "ENSG00000128604", gene_name = "IRF5", strand = 1))

#3'UTR QTL
#Construct metadata df for wiggleplotr
samples_in_dir = data_frame(bigWig = dir("/Volumes/Ajamasin/bigwig/RNA//")) %>% 
  tidyr::separate(bigWig, c("sample_id", "strand","suffix"), sep = "\\.")
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, "rs10954213", vcf_file$genotypes, -1)
track_data = dplyr::semi_join(track_data, samples_in_dir, by = "sample_id")

#Filter by condtion
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive","IFNg"))

#Make a coverage plot of the ATAC data
down_exons = c(union_exons, alt_events$ENSG00000128604.clique_1$downstream)
down_cdss = c(union_cds, alt_events$ENSG00000128604.clique_1$downstream)
IRF5_utr_plot = plotCoverage(exons = down_exons, 
                            cdss = down_cdss, 
                            track_data = filtered_tracks, rescale_introns = TRUE, 
                            transcript_annotations = union_annotations, fill_palette = getGenotypePalette(), 
                            plot_fraction = 0.2, heights = c(0.55,0.45), 
                            return_subplots_list = FALSE, line_only = TRUE)
ggsave("figures/supplementary/IRF5_UTR_qtl.pdf", plot = IRF5_utr_plot, width = 9, height = 6)


# Alternative promoter QTL
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, "rs3778754", vcf_file$genotypes, -1)
track_data = dplyr::semi_join(track_data, samples_in_dir, by = "sample_id")

#Filter by condtion
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive","IFNg"))

#Make a coverage plot of the ATAC data
up_exons = c(union_exons, alt_events$ENSG00000128604.clique_1$upstream)
up_cdss = c(union_cds, alt_events$ENSG00000128604.clique_1$upstream)
IRF5_promoter_plot = plotCoverage(exons = up_exons, 
                             cdss = up_cdss, 
                             track_data = filtered_tracks, rescale_introns = TRUE, 
                             transcript_annotations = union_annotations, fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.55,0.45), 
                             return_subplots_list = FALSE, line_only = TRUE)
ggsave("figures/supplementary/IRF5_promoter_qtl.pdf", plot = IRF5_promoter_plot, width = 9, height = 6.5)

# RASQUAL QTL
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, "rs199508964", vcf_file$genotypes, -1)
track_data = dplyr::semi_join(track_data, samples_in_dir, by = "sample_id")

#Filter by condtion
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive","IFNg"))

#Make a coverage plot of the ATAC data
up_exons = c(union_exons, as.list(gene_data_ext$exons)[c("ENST00000619830","ENST00000402030","ENST00000477535")])
up_cdss = c(union_cds, as.list(gene_data_ext$cds)[c("ENST00000619830","ENST00000402030","ENST00000477535")])
IRF5_rasqual_plot = plotCoverage(exons = up_exons, 
                                  cdss = up_cdss, 
                                  track_data = filtered_tracks, rescale_introns = TRUE, 
                                  transcript_annotations = union_annotations, fill_palette = getGenotypePalette(), 
                                  plot_fraction = 0.2, heights = c(0.55,0.45), 
                                  return_subplots_list = FALSE, line_only = TRUE)
ggsave("figures/supplementary/IRF5_rasqual_qtl.pdf", plot = IRF5_rasqual_plot, width = 9, height = 6.5)



#Make QTL plots for the IRF5 variants
se_ensembl = readRDS("results/SL1344/combined_reviseAnnotations_transcript_quants.rds")
ratio_matrix = assays(se_ensembl)$tpm_ratios
normalised_matrix = ratio_matrix %>% replaceNAsWithRowMeans() %>% quantileNormaliseRows()
sample_meta = colData(se_ensembl) %>% tbl_df2() %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344")))
gene_meta = rowData(se_ensembl) %>% tbl_df2() %>%
  dplyr::mutate(gene_id = transcript_id)

#Import Salmon QTLs
salmon_qtl_hits = readRDS("results/SL1344/salmon/salmon_qtl_hits.rds")
salmon_qtl_df = purrr::map_df(salmon_qtl_hits, identity, .id = "condition_name")
dplyr::filter(salmon_qtl_df, gene_name == "IRF5")

#Make plots on raw data
qtl_df_utr = constructQtlPlotDataFrame("ENSG00000128604.clique_1.downstream.ENST00000489702", "rs10954213", 
                                   ratio_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::filter(norm_exp < 0.75) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  plotQtlRow(ylabel = "Relative expression")
ggsave("figures/supplementary/IRF5_UTR_boxplot.pdf", plot = qtl_df_utr, width = 4.5, height = 3)
qtl_df_promoter = constructQtlPlotDataFrame("ENSG00000128604.clique_1.upstream.ENST00000249375", "rs3778754", 
                                   ratio_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  plotQtlRow()
ggsave("figures/supplementary/IRF5_promoter_boxplot.pdf", plot = qtl_df_promoter, width = 4.5, height = 3)

#Fetch SNPs around the IRF5 gene
fastqtl_meta = dplyr::transmute(gene_meta, gene_id, chr, start = transcript_start, end = transcript_end)
region_coords = c(128935000,128952000)
region_coords = c(128900000,128990000)

utr_pvalues = fastqtlTabixFetchGenesQuick("ENSG00000128604.clique_1.downstream.ENST00000489702", 
                                          "/Volumes/Ajamasin/databases/SL1344/salmon/naive_100kb_pvalues.sorted.txt.gz", 
                                          gene_metadata = fastqtl_meta, cis_window = 1e5)[[1]] %>%
  dplyr::mutate(track_id = "naive") %>% 
  dplyr::select(gene_id, pos, p_nominal, track_id)
promoter_pvalues = fastqtlTabixFetchGenesQuick("ENSG00000128604.clique_1.upstream.ENST00000249375", 
                                          "/Volumes/Ajamasin/databases/SL1344/salmon/naive_100kb_pvalues.sorted.txt.gz", 
                                          gene_metadata = fastqtl_meta, cis_window = 1e5)[[1]] %>%
  dplyr::mutate(track_id = "naive") %>% 
  dplyr::select(gene_id, pos, p_nominal, track_id)

#Fetch rasqual SNPs around the lead variant
rasqual_pvalues = tabixFetchGenesQuick("ENSG00000128604", qtlResults()$rna_rasqual$naive, 
                            gene_metadata = combined_expression_data$gene_metadata, cis_window = 1e5)[[1]] %>%
  dplyr::mutate(track_id = "naive") %>% 
  dplyr::select(gene_id, pos, p_nominal, track_id)

utr_manhatten = makeManhattanPlot(utr_pvalues, region_coords)
promoter_manhattan = makeManhattanPlot(promoter_pvalues, region_coords)
rasqual_manhattan = makeManhattanPlot(rasqual_pvalues, region_coords)

#Calculate the union of exons
union_exons = list(ENST00000467002 = listUnion(gene_data_ext$exons))
union_cds = list(ENST00000467002 = listUnion(gene_data_ext$cdss))
union_transcripts = wiggleplotr::plotTranscripts(union_exons, union_cds, plotting_annotations, rescale_introns = FALSE, region_coords = region_coords)

promoter_qtl_plot = cowplot::plot_grid(promoter_manhattan, union_transcripts, 
    align = "v", rel_heights = c(0.6,0.4), ncol = 1)
utr_qtl_plot = cowplot::plot_grid(utr_manhatten, union_transcripts, 
                   align = "v", rel_heights = c(0.6,0.4), ncol = 1)
rasqual_qtl_plot = cowplot::plot_grid(rasqual_manhattan, union_transcripts, 
                                  align = "v", rel_heights = c(0.6,0.4), ncol = 1)
ggsave("figures/supplementary/IRF5_promoter_manhattan.pdf", plot = promoter_qtl_plot, width = 9, height = 3)
ggsave("figures/supplementary/IRF5_UTR_manhattan.pdf", plot = utr_qtl_plot, width = 9, height = 3)
ggsave("figures/supplementary/IRF5_rasqual_manhattan.pdf", plot = rasqual_qtl_plot, width = 9, height = 3)


#Look at the LeafCutter QTLs as well
leafcutter_min_pvalues = readRDS("results/SL1344/leafcutter/leafcutter_cluster_min_pvalues.rds")
leafcutter_gene_pvalues = purrr::map(leafcutter_min_pvalues, ~dplyr::arrange(., ensembl_gene_id, p_nominal) %>%
                                       dplyr::filter(!is.na(ensembl_gene_id)) %>%
                                       dplyr::group_by(ensembl_gene_id) %>% 
                                       dplyr::filter(row_number() == 1) %>%
                                       dplyr::ungroup())
purrr::map_df(leafcutter_gene_pvalues, identity, .id = "condition_name") %>% dplyr::filter(ensembl_gene_id == "ENSG00000128604")

prop_list = readRDS("results/SL1344/combined_proportions.row_quantile.rds")
cluster_meta = dplyr::select(prop_list$gene_metadata, gene_id, cluster_id, cluster_size)                                     
                                     

wiggleplotr::plotTranscripts(gene_data_ext$exons, gene_data_ext$cdss, plotting_annotations, rescale_introns = FALSE, region_coords = c(128946056,128948278))


