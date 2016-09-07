library("GenomicFeatures")
library("dplyr")
library("devtools")
load_all("../reviseAnnotations/")
load_all("../wiggleplotr/")
load_all("../seqUtils/")

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

#Extract gene data from annotations
gene_data = reviseAnnotations::extractGeneData("ENSG00000003056", filtered_metadata, exons, cdss)

#Extend truncated transcripts until the longest transcript
gene_extended_tx = reviseAnnotations::extendTranscriptsPerGene(gene_data$metadata, gene_data$exons, gene_data$cdss)
gene_data_ext = reviseAnnotations::replaceExtendedTranscripts(gene_data, gene_extended_tx)

a = constructAlternativeEvents(gene_data_ext$exons, "ENSG00000103495")

#Make some plots
plotting_annotations = select(filtered_metadata, ensembl_transcript_id, ensembl_gene_id, external_gene_name, strand) %>% 
  dplyr::rename(transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)

wiggleplotr::plotTranscripts(gene_data$exons, gene_data$cdss, plotting_annotations, rescale_introns = TRUE)
wiggleplotr::plotTranscripts(gene_data_ext$exons, gene_data_ext$cdss, plotting_annotations, rescale_introns = TRUE)
wiggleplotr::plotTranscripts(a$ENSG00000103495.clique_1$contained, a$ENSG00000103495.clique_1$contained, plotting_annotations, rescale_introns = TRUE)

exons1 = findTranscriptCliques(gene_data_ext$exons)[[2]]
wiggleplotr::plotTranscripts(gene_data_ext$exons[exons1], gene_data_ext$cdss[exons1], plotting_annotations, rescale_introns = TRUE)

exons1 = findTranscriptCliques(gene_data$exons)[[2]]
wiggleplotr::plotTranscripts(gene_data$exons[exons1], gene_data$cdss[exons1], plotting_annotations, rescale_introns = TRUE)
