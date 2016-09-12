library("GenomicFeatures")
library("dplyr")
library("devtools")
load_all("../reviseAnnotations/")
load_all("../wiggleplotr/")
load_all("../seqUtils/")
library("ggplot2")

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



#Make coverage plots
#Construct metadata df for wiggleplotr
samples_in_dir = data_frame(bigWig = dir("/Volumes/Ajamasin/bigwig/RNA//")) %>% 
  tidyr::separate(bigWig, c("sample_id", "strand","suffix"), sep = "\\.")
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str1_df, "rs10954213", vcf_file$genotypes, 1)
track_data = dplyr::semi_join(track_data, samples_in_dir, by = "sample_id")

#Filter by condtion
filtered_tracks = dplyr::filter(track_data, track_id %in% c("IFNg"))

#Make a coverage plot of the ATAC data
selected_transcripts = c("ENST00000473745","ENST00000249375","ENST00000489702","ENST00000357234")
IRF5_utr_plot = plotCoverage(exons = gene_data_ext$exons[selected_transcripts], 
                            cdss = gene_data_ext$cds[selected_transcripts], 
                            track_data = filtered_tracks, rescale_introns = TRUE, 
                            transcript_annotations = plotting_annotations, fill_palette = getGenotypePalette(), 
                            plot_fraction = 0.2, heights = c(0.7,0.3), 
                            return_subplots_list = FALSE)
ggsave("figures/supplementary/IRF5_UTR_qtl.pdf", plot = IRF5_utr_plot, width = 7, height = 6)



track_data = wiggleplotrGenotypeColourGroup(rna_meta_str1_df, "rs3778754", vcf_file$genotypes, -1)
track_data = dplyr::semi_join(track_data, samples_in_dir, by = "sample_id")

#Filter by condtion
filtered_tracks = dplyr::filter(track_data, track_id %in% c("IFNg"))

#Make a coverage plot of the ATAC data
selected_transcripts = c("ENST00000473745","ENST00000249375","ENST00000489702","ENST00000357234")
IRF5_promoter_plot = plotCoverage(exons = gene_data_ext$exons[selected_transcripts], 
                             cdss = gene_data_ext$cds[selected_transcripts], 
                             track_data = filtered_tracks, rescale_introns = TRUE, 
                             transcript_annotations = plotting_annotations, fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.7,0.3), 
                             return_subplots_list = FALSE)
ggsave("figures/supplementary/IRF5_promoter_qtl.pdf", plot = IRF5_promoter_plot, width = 7, height = 6)









