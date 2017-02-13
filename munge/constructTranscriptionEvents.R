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
gene_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.filtered.rds")  
txdb = AnnotationDbi::loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
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

#Construct alternative events and remove failed genes
safe_construct = purrr::safely(constructAlternativeEventsWrapper)
alt_events = purrr::map(gene_ids_list, ~safe_construct(., filtered_metadata, exons, cdss)$result)
failed_genes = purrr::map_lgl(alt_events, is.null)
alt_events = alt_events[!failed_genes] #Remove failed genes

#Flatten
alt_events = purrr::flatten(alt_events) %>% flattenAlternativeEvents()

#Construct event metadata
event_metadata = reviseAnnotations::constructEventMetadata(names(alt_events))

#Construct transcript annotations
transcript_annotations = reviseAnnotations::transcriptsToAnnotations(alt_events, event_metadata)

#Make a list of failed genes
failed_names = names(which(failed_genes))

#Save output from each batch
if(!is.null(batch_id)){
  output_file = file.path("results/reviseAnnotations", paste0("reviseAnnotations_batch_",batch_id, ".gff3"))
  error_file = file.path("results/reviseAnnotations", paste0("failed_genes_batch_",batch_id, ".txt"))
  rtracklayer::export.gff3(transcript_annotations, output_file)
  write.table(failed_names, error_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
