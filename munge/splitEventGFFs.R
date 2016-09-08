library("rtracklayer")
library("dplyr")
library("devtools")
load_all("../reviseAnnotations/")
load_all("../seqUtils/")

#Imort the full GFF file
gff = rtracklayer::import.gff3("results/reviseAnnotations/reviseAnnotations.gff3.gz")

#Construct metadata
tx_records = dplyr::filter(as.data.frame(gff) %>% tbl_df(), type == "mRNA") %>%
  dplyr::transmute(transcript_id = ID, transcript_start = start, transcript_end = end, 
                   transcript_length = width, strand = as.character(strand), chr = as.character(seqnames))
event_metadata = reviseAnnotations::constructEventMetadata(tx_records$transcript_id)
tx_metadata = dplyr::left_join(tx_records, event_metadata, by = "transcript_id")
saveRDS(tx_metadata, "results/reviseAnnotations/reviseAnnotations.transcript_metadata.rds")

#split events by position
upstream_gff = gff[grep("upstream", gff$gene_id),]
downstream_gff = gff[grep("downstream", gff$gene_id),]
contained_gff = gff[grep("contained", gff$gene_id),]

#Export different events
rtracklayer::export.gff3(upstream_gff, "results/reviseAnnotations/reviseAnnotations.upstream.gff3")
rtracklayer::export.gff3(downstream_gff, "results/reviseAnnotations/reviseAnnotations.downstream.gff3")
rtracklayer::export.gff3(contained_gff, "results/reviseAnnotations/reviseAnnotations.contained.gff3")
