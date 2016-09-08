library("DRIMSeq")
library("SummarizedExperiment")
library("devtools")
load_all("../seqUtils/")

#Import event dataset
se_ensembl = readRDS("results/SL1344/combined_reviseAnnotations_transcript_quants.rds")
event_metadata = rowData(se_ensembl)
unique_genes = unique(event_metadata$ensembl_gene_id)

####### Get batch string from stdin ######
f <- file("stdin")
open(f)
batch_string = readLines(f) %>% as.numeric()
close(f)
####### END #######

#Construct splicing dataset
splicing_data <- DRIMSeq::dmDSdata(counts = assays(se_ensembl)$counts, 
                                   gene_id = rowData(se_ensembl)$gene_id,
                                   feature_id = rowData(se_ensembl)$transcript_id, 
                                   sample_id = se_ensembl$sample_id, 
                                   group = se_ensembl$condition_name)

#Filter by counts
splicing_data_filtered <- dmFilter(splicing_data, min_samps_gene_expr = 84, min_samps_feature_expr = 84,
                          min_samps_feature_prop = 0)

#Select genes from a single batch
batch_genes = constructIdBatches(batch_string, unique_genes)
selected_genes = dplyr::filter(rowData(se_ensembl) %>% as.data.frame(), ensembl_gene_id %in% batch_genes)$gene_id
splicing_data_batch <- splicing_data_filtered[names(splicing_data_filtered) %in% selected_genes, ]

#Run DRIMseq
sd <- dmDispersion(splicing_data_batch, verbose = 1, BPPARAM = BiocParallel::SerialParam())
sd <- dmFit(sd, BPPARAM = BiocParallel::SerialParam())

#Construct output data frame
prop_df = proportions(sd)

#Perform all three tests
SL1344_dte <- dmTest(sd, verbose = 1, BPPARAM = BiocParallel::SerialParam(), compared_groups = c(1,3)) %>%
  results() %>% dplyr::transmute(gene_id, SL1344_lr = lr, SL1344_pvalue = pvalue)
IFNg_dte <- dmTest(sd, verbose = 1, BPPARAM = BiocParallel::SerialParam(), compared_groups = c(1,2)) %>%
  results() %>% dplyr::transmute(gene_id, IFNg_lr = lr, IFNg_pvalue = pvalue)
IFNg_SL1344_dte <- dmTest(sd, verbose = 1, BPPARAM = BiocParallel::SerialParam(), compared_groups = c(1,3)) %>%
  results() %>% dplyr::transmute(gene_id, IFNg_SL1344_lr = lr, IFNg_SL1344_pvalue = pvalue)
all_dte <- dmTest(sd, verbose = 1, BPPARAM = BiocParallel::SerialParam()) %>%
  results() %>% dplyr::transmute(gene_id, all_lr = lr, all_pvalue = pvalue)

#Merge all results
results_df = dplyr::left_join(prop_df, SL1344_dte, by = "gene_id") %>% 
  dplyr::left_join(IFNg_dte, by = "gene_id") %>%
  dplyr::left_join(IFNg_SL1344_dte, by = "gene_id") %>%
  dplyr::left_join(all_dte, by = "gene_id")

#Save output from each batch
if(!is.null(batch_id)){
  output_file = file.path("results/reviseAnnotations", paste0("differential_events_batch",gsub(" ","_",batch_string), ".gff3"))
  write.table(results_df, output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
