library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import sample metadata
gene_data = readRDS("results/SL1344/combined_expression_data.rds")
sample_metadata = gene_data$sample_metadata %>% as.data.frame()
rownames(sample_metadata) = sample_metadata$sample_id

#Import transcript metadata
tx_meta = readRDS("../../annotations/GRCh38/genes/Ensembl_85/Homo_sapiens.GRCh38.85.transcript_data.rds") %>%
  dplyr::rename(gene_id = ensembl_gene_id, gene_name = external_gene_name) %>%
  dplyr::mutate(transcript_id = paste(ensembl_transcript_id, transcript_version, sep = "."))

#Construct SummarizedExeperiment object
se_ensembl = salmonSummarizedExperiment(sample_metadata, tx_meta, "results/SL1344/salmon_quants/", sub_dir = FALSE)
saveRDS(se_ensembl,"results/SL1344/combined_ensembl85_transcript_quants.rds")