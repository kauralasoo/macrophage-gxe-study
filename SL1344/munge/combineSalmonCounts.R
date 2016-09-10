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
  dplyr::mutate(transcript_id = paste(ensembl_transcript_id, transcript_version, sep = "."), ensembl_gene_id = gene_id)

#Construct SummarizedExeperiment object
se_ensembl = salmonSummarizedExperiment(sample_metadata, tx_meta, "results/SL1344/salmon_quants/", 
                                        sub_dir = FALSE, counts_suffix = ".ensembl85.quant.sf.gz", skip = 1)
saveRDS(se_ensembl,"results/SL1344/combined_ensembl85_transcript_quants.rds")

#Import Salmon event quantifications from reviseAnnotations
transcript_meta = readRDS("results/reviseAnnotations/reviseAnnotations.transcript_metadata.rds") %>%
  dplyr::left_join(dplyr::select(gene_data$gene_metadata, gene_id, gene_name), by = c("ensembl_gene_id" = "gene_id"))
se_reviseAnnotations = salmonSummarizedExperiment(sample_metadata, transcript_meta, "results/SL1344/salmon_reviseAnnotations_quants/", 
                                                  counts_suffix = ".reviseAnnotations.quant.sf.gz", sub_dir = FALSE)
saveRDS(se_reviseAnnotations,"results/SL1344/combined_reviseAnnotations_transcript_quants.rds")
