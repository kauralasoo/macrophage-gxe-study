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
gene_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.rds")  
txdb = AnnotationDbi::loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)

#Filter gene annotations
filtered_metadata = reviseAnnotations::filterTranscriptMetadata(gene_metadata)

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

#Separate the two groups
grp1_events = dplyr::filter(event_metadata, grp_id == "grp_1")
grp2_events = dplyr::filter(event_metadata, grp_id == "grp_2")

#Construct transcript annotations
grp1_annotations = reviseAnnotations::transcriptsToAnnotations(alt_events[grp1_events$transcript_id], grp1_events)
grp2_annotations = reviseAnnotations::transcriptsToAnnotations(alt_events[grp2_events$transcript_id], grp2_events)

#Make a list of failed genes
failed_names = names(which(failed_genes))

#Save output from each batch
if(!is.null(batch_id)){
  grp1_file = file.path("results/reviseAnnotations", paste0("reviseAnnotations.grp_1.batch_",batch_id, ".gff3"))
  grp2_file = file.path("results/reviseAnnotations", paste0("reviseAnnotations.grp_2.batch_",batch_id, ".gff3"))
  error_file = file.path("results/reviseAnnotations", paste0("failed_genes.batch_",batch_id, ".txt"))
  rtracklayer::export.gff3(grp1_annotations, grp1_file)
  rtracklayer::export.gff3(grp2_annotations, grp2_file)
  write.table(failed_names, error_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
