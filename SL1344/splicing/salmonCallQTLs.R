library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("Rsamtools")
library("purrr")
library("SummarizedExperiment")

#Import event dataset
se_ensembl = readRDS("results/SL1344/combined_reviseAnnotations_transcript_quants.rds")
transcript_metadata = rowData(se_ensembl) %>% tbl_df2() %>% 
  dplyr::select(transcript_id, gene_id, event_type, ensembl_gene_id, gene_name)

#Identify expressed genes
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
mean_expression = calculateMean(combined_expression_data$tpm, as.data.frame(combined_expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_expression, 1, max) > .5))

#Find minimal p-values from fastQTL results
naive_fqtl = importFastQTLTable("results/SL1344/salmon/fastqtl_output/naive_100kb_permuted.txt.gz") %>%
  dplyr::rename(transcript_id = gene_id)
ifng_fqtl = importFastQTLTable("results/SL1344/salmon/fastqtl_output/IFNg_100kb_permuted.txt.gz") %>%
  dplyr::rename(transcript_id = gene_id)
sl1344_fqtl = importFastQTLTable("results/SL1344/salmon/fastqtl_output/SL1344_100kb_permuted.txt.gz") %>%
  dplyr::rename(transcript_id = gene_id)
ifng_sl1344_fqtl = importFastQTLTable("results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_permuted.txt.gz") %>%
  dplyr::rename(transcript_id = gene_id)

#Construct a list
fastqtl_pvalue_list = list(naive = naive_fqtl,
                           IFNg = ifng_fqtl,
                           SL1344 = sl1344_fqtl, 
                           IFNg_SL1344 = ifng_sl1344_fqtl)

#Remove duplicate transcript entries (do not know why they appeared)
fastqtl_pvalue_unique = purrr::map(fastqtl_pvalue_list, ~dplyr::group_by(.,transcript_id) %>% 
                                     dplyr::filter(row_number() == 1) %>% ungroup())

#Add transcript metadata to QTLs
fastqtl_pvalue_meta = purrr::map(fastqtl_pvalue_unique, ~dplyr::left_join(.,transcript_metadata, by = "transcript_id"))
saveRDS(fastqtl_pvalue_meta, "results/SL1344/salmon/salmon_min_pvalues_meta.rds")

#Keep only expressed genes
fastqtl_pvalue_expressed = purrr::map(fastqtl_pvalue_meta, ~dplyr::filter(.,ensembl_gene_id %in% expressed_genes))

#Identify transcription QTLs
salmon_qtl_hits = purrr::map(fastqtl_pvalue_expressed, ~dplyr::group_by(., ensembl_gene_id, event_type) %>% 
                               dplyr::mutate(n_transcripts = length(transcript_id)) %>% 
                               dplyr::arrange(ensembl_gene_id, p_beta) %>% dplyr::filter(row_number() == 1) %>% 
                               dplyr::ungroup() %>% dplyr::mutate(p_bonferroni = p_beta * n_transcripts) %>% 
                               dplyr::mutate(p_bonferroni = pmin(p_bonferroni, 1)) %>% 
                               dplyr::mutate(p_fdr = p.adjust(p_bonferroni, method = "fdr")) %>% 
                               dplyr::filter(p_fdr < 0.1))
saveRDS(salmon_qtl_hits, "results/SL1344/salmon/salmon_qtl_hits.rds")
salmon_qtl_df = purrr::map_df(salmon_qtl_hits, identity, .id = "condition_name")






