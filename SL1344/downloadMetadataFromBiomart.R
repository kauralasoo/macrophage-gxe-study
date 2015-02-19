library("biomaRt")
library("dplyr")

ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", host = "dec2014.archive.ensembl.org")
ensembl_dataset = useDataset("hsapiens_gene_ensembl",mart=ensembl_mart)

#Download data from biomart
selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name" ,"gene_biotype", 
                        "transcript_biotype","chromosome_name","strand")
data = getBM(attributes = selected_attributes, mart = ensembl_dataset)

#Rename some columns
data = dplyr::rename(data, transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)
saveRDS(data, "annotations/biomart_transcripts.rds")

#Create a filtered version of the annotations
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")

#Filter the data set
filtered_data = dplyr::filter(data, chromosome_name %in% valid_chromosomes, 
                              gene_biotype %in% valid_gene_biotypes,
                              transcript_biotype %in% valid_gene_biotypes)
saveRDS(filtered_data, "annotations/biomart_transcripts.filtered.rds")

