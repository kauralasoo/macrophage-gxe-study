library("dplyr")
library("devtools")
library("GenomicFeatures")
library("wiggleplotr")
load_all("../seqUtils/")
load_all("../reviseAnnotations/")



#Import transcript annotations
gene_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.rds")  
txdb = AnnotationDbi::loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)

#Filter gene annotations
filtered_metadata = reviseAnnotations::filterTranscriptMetadata(gene_metadata)

#Extract gene data from annotations
gene_data = reviseAnnotations::extractGeneData("ENSG00000128604", filtered_metadata, exons, cdss)

#Extend truncated transcripts until the longest transcript
gene_extended_tx = reviseAnnotations::extendTranscriptsPerGene(gene_data$metadata, gene_data$exons, gene_data$cdss)
gene_data_ext = reviseAnnotations::replaceExtendedTranscripts(gene_data, gene_extended_tx)

alt_events = constructAlternativeEvents(gene_data_ext$exons, "ENSG00000128604")

#Make some plots
plotting_annotations = dplyr::select(filtered_metadata, ensembl_transcript_id, ensembl_gene_id, external_gene_name, strand) %>% 
  dplyr::rename(transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)

wiggleplotr::plotTranscripts(gene_data$exons, gene_data$cdss, plotting_annotations, rescale_introns = TRUE)
wiggleplotr::plotTranscripts(gene_data_ext$exons, gene_data_ext$cdss, plotting_annotations, rescale_introns = TRUE)
wiggleplotr::plotTranscripts(alt_events[[1]]$upstream, alt_events[[1]]$upstream, plotting_annotations, rescale_introns = TRUE)
wiggleplotr::plotTranscripts(alt_events[[1]]$downstream, alt_events[[1]]$downstream, plotting_annotations, rescale_introns = TRUE)
wiggleplotr::plotTranscripts(alt_events[[1]]$contained, alt_events[[1]]$contained, plotting_annotations, rescale_introns = TRUE)

wiggleplotr::plotTranscripts(alt_events[[2]]$upstream, alt_events[[2]]$upstream, plotting_annotations, rescale_introns = TRUE)
wiggleplotr::plotTranscripts(alt_events[[2]]$downstream, alt_events[[2]]$downstream, plotting_annotations, rescale_introns = TRUE)
wiggleplotr::plotTranscripts(alt_events[[2]]$contained, alt_events[[2]]$contained, plotting_annotations, rescale_introns = TRUE)

