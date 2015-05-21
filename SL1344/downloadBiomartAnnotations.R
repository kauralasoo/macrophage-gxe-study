library("GenomicFeatures")
library("biomaRt")
library("dplyr")

#Make TranscriptDb object from biomart
#txdb78 = makeTranscriptDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="dec2014.archive.ensembl.org")
#saveDb(txdb78, "annotations/TranscriptDb_GRCh38_78.db")
#exons = exonsBy(txdb78, by = "tx", use.names = TRUE)

#Downlaod transcript metadata from Ensembl
ensembl = useMart("ENSEMBL_MART_ENSEMBL", host="may2015.archive.ensembl.org")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)

#Define attributes to be downloaded from biomart
biomart_attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype","status", 
                       "chromosome_name","strand", "transcript_start","transcript_end","ensembl_transcript_id", 
                       "transcript_status", "transcript_tsl",
                       "transcript_gencode_basic", "transcript_appris_pi","external_transcript_name", 
                       "transcript_length", "transcript_biotype", "ccds")
refseq_attributes = c("ensembl_gene_id","ensembl_transcript_id","refseq_mrna")

#Download data for sample genes
genes = c("ENSG00000111912", "ENSG00000266094","ENSG00000068028")
data = getBM(attributes = biomart_attributes, 
             filters = c("ensembl_gene_id"), 
             values = genes, 
             mart = ensembl)
data

data1 = getBM(attributes = refseq_attributes, 
              filters = c("ensembl_gene_id"), 
              values = genes, 
              mart = ensembl)
data1

#Download data for all transcripts
transcript_data = getBM(attributes = biomart_attributes, mart = ensembl)
refseq_data = getBM(attributes = refseq_attributes, mart = ensembl)
saveRDS(transcript_data, "annotations/Homo_sapiens.GRCh38.78.transcript_data.rds")
saveRDS(refseq_data, "annotations/Homo_sapiens.GRCh38.78.refseq_data.rds")



