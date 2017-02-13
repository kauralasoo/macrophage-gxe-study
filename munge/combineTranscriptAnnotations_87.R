library("dplyr")
library("tidyr")
library("devtools")
load_all("../seqUtils/")
load_all("../reviseAnnotations/")

#Read different datasets from disk
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.transcript_data.rds"))
refseq_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.refseq_data.rds"))
transcript_tags = read.table("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.transcript_tags.txt", stringsAsFactors = FALSE)

#Extract transcript tags
colnames(transcript_tags) = c("ensembl_transcript_id", "tags")
tag_list = strsplit(transcript_tags$tags, ";")
names(tag_list) = transcript_tags$ensembl_transcript_id

#Build a data.frame of transcript tags
tags_df = transmute(transcript_data, ensembl_transcript_id,  CCDS = 0, mRNA_start_NF = 0, mRNA_end_NF = 0, cds_start_NF = 0, cds_end_NF = 0, seleno = 0)
tags_df = data.frame(tags_df)
rownames(tags_df) = tags_df$ensembl_transcript_id

#Put tags into the data frame
#NOTE: This will take a while...
tx_ids = names(tag_list)
for (tx_id in tx_ids){
  print(tx_id)
  tags_df[tx_id, tag_list[[tx_id]] ] = 1
}
saveRDS(tbl_df(tags_df), "../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.tags_dataframe.rds")

#Add transcript tags into the data.frame
tags_df = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.tags_dataframe.rds")
tags_df = dplyr::select(tags_df, ensembl_transcript_id, mRNA_start_NF:cds_end_NF)

#Convert TSL to integer
trimmed_tsl = dplyr::select(transcript_data, ensembl_transcript_id, transcript_tsl) %>% 
  tidyr::separate(transcript_tsl, into = c("transcript_tsl", "version"), sep = "\\s\\(", extra = "drop") %>% 
  dplyr::select(ensembl_transcript_id, transcript_tsl) %>%
  dplyr::mutate(transcript_tsl = ifelse(transcript_tsl == "", "tslNA", transcript_tsl)) %>% #Convert empty TSLs into NAs 
  tidyr::separate(transcript_tsl, into = c("none","tsl"), sep ="tsl", extra = "drop", convert = TRUE) %>% 
  dplyr::select(ensembl_transcript_id, tsl)

#Add refseq ids 
refseq_ids = dplyr::mutate(refseq_data, refseq_id_count = ifelse(refseq_mrna == "", 0, 1)) %>% 
  dplyr::group_by(ensembl_transcript_id) %>% 
  dplyr::summarise(refseq_id_count = sum(refseq_id_count), refseq_ids = paste(refseq_mrna, collapse = ";"))

#Put all of the data together
compiled_data = dplyr::left_join(transcript_data, tags_df, by = "ensembl_transcript_id") %>% #Add transcript flags
  dplyr::mutate(mRNA_start_end_NF = mRNA_start_NF + mRNA_end_NF) %>%
  dplyr::mutate(cds_start_end_NF = cds_start_NF + cds_end_NF) %>%
  dplyr::left_join(trimmed_tsl, by = "ensembl_transcript_id") %>% #Add TSL 
  dplyr::left_join(refseq_ids, by = "ensembl_transcript_id") %>% #Add RefSeq ids
  dplyr::mutate(is_ccds = ifelse(ccds == "", 0, 1)) %>%
  dplyr::mutate(is_gencode_basic = ifelse(transcript_gencode_basic == "",0,1))

#Count CCDS and GENCODE basic annotations per gene
compiled_data = dplyr::group_by(compiled_data, ensembl_gene_id) %>% 
  dplyr::summarize(n_ccds = sum(is_ccds), n_gencode_basic = sum(is_gencode_basic)) %>% 
  dplyr::left_join(compiled_data, by = "ensembl_gene_id")
saveRDS(compiled_data, "../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.rds")  

#Filter annotations
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")

#Remove retained introns from the annotations
filtered_data = dplyr::filter(compiled_data, chromosome_name %in% valid_chromosomes, 
                              gene_biotype %in% valid_gene_biotypes,
                              transcript_biotype %in% valid_gene_biotypes) %>%
  #Flag transcripts that are potenitally good references
  dplyr::mutate(is_good_reference = ifelse(cds_start_end_NF == 0, is_gencode_basic, 0)) 
saveRDS(filtered_data, "../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.filtered.rds")  

#For each gene mark the transcripts with longest starts and ends
marked_data = markLongestTranscripts(filtered_data)
saveRDS(marked_data, "../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.longest_marked.rds")  

#Are multiple Ids for the same gene a problem
multiple_ids_per_name = dplyr::select(filtered_data, ensembl_gene_id, external_gene_name, gene_biotype) %>% 
  unique() %>% 
  group_by(external_gene_name) %>% 
  summarize(id_count = length(ensembl_gene_id), gene_biotype = gene_biotype[1]) %>% 
  arrange(desc(id_count)) %>%
  filter(id_count > 1)







