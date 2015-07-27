library("plyr")
library("dplyr")
library("tidyr")
library("devtools")
library("GenomicFeatures")
load_all("macrophage-gxe-study/seqUtils/")

#Load transcript annotations from disk
txdb79 = loadDb("../../annotations/GRCh38/genes/TranscriptDb_GRCh38_79.db")
exons = exonsBy(txdb79, by = "tx", use.names = TRUE)

#Load Gencode BASIC transcripts
gencode_basic_transcript = read.table("../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.txt", header = TRUE, stringsAsFactors = FALSE) %>% 
  tbl_df() %>%
  dplyr::rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id)

#Convert transcript annotations to 
gene_exon_start_end = extractExonsStartEnd(exons, gencode_basic_transcript)
write.table(gene_exon_start_end, "annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
