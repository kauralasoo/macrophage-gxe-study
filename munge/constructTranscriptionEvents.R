library("GenomicFeatures")
library("dplyr")
library("devtools")
load_all("../reviseAnnotations/")
load_all("../seqUtils/")

batch_id = NULL

####### Get batch id from disk ######
f <- file("stdin")
open(f)
batch_id = readLines(f) %>% as.numeric()
close(f)
####### END #######

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

#### Split genes into batches ####
gene_ids = unique(filtered_metadata$ensembl_gene_id)
batches = splitIntoBatches(length(gene_ids), 50)
if(!is.null(batch_id)){
  selection = batches == batch_id
  gene_ids = gene_ids[selection]
}

#Construct events
gene_ids_list = seqUtils::idVectorToList(gene_ids)
alt_events = purrr::map(gene_ids_list, ~constructAlternativeEventsWrapper(., filtered_metadata, exons, cdss)) %>% 
  purrr::flatten() %>%
  flattenAlternativeEvents()

#Construct event metadata
event_metadata = data.frame(transcript_id = names(alt_events)) %>% 
  tidyr::separate(transcript_id, c('ensembl_gene_id', 'clique_id', 'event_type','ensembl_transcript_id'), 
                  sep = "\\.", remove = F) %>%
  dplyr::mutate(gene_id = paste(ensembl_gene_id, clique_id, event_type, sep = ".")) %>%
  dplyr::mutate(transcript_id = as.character(transcript_id))

#Construct transcript annotations
transcript_annotations = transcriptsToAnnotations(alt_events, event_metadata)

#Save output from each batch
if(!is.null(batch_id)){
  output_file = file.path("results/reviseAnnotations", paste0("/reviseAnnotations_batch_",batch_id, ".gff3"))
  rtracklayer::export.gff3(transcript_annotations, output_file)
}
