library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("Rsamtools")
library("purrr")

#Import leafCutter data
leafcutter_se = readRDS("results/SL1344/leafCutter_summarized_experiment.rds")
event_metadata = rowData(leafcutter_se) %>% tbl_df2()

#Import gene metadata
combined_expression_data_filtered = readRDS("results/SL1344/combined_expression_data.rds")
gene_ranges = dplyr::transmute(combined_expression_data_filtered$gene_metadata, gene_id, seqnames = chr, start, end, strand = "+") %>% 
  dataFrameToGRanges()

#Link leafCutter clusters to genes
cluster_ranges = event_metadata %>% 
  dplyr::transmute(cluster_id, seqnames = chr, start = cluster_start, end = cluster_end, strand = "+") %>% 
  unique() %>% 
  dataFrameToGRanges()
olaps = findOverlaps(cluster_ranges, gene_ranges)

cluster_gene_map = data_frame(cluster_id = cluster_ranges[queryHits(olaps),]$cluster_id, ensembl_gene_id = gene_ranges[subjectHits(olaps),]$gene_id) %>% 
  dplyr::group_by(cluster_id) %>% 
  dplyr::mutate(overlap_count = length(ensembl_gene_id)) %>% 
  dplyr::filter(overlap_count == 1) %>% 
  dplyr::select(-overlap_count)
event_metadata_names = dplyr::left_join(event_metadata, cluster_gene_map, by = "cluster_id")

#Find minimal p-values from fastQTL results
naive_fqtl = importFastQTLTable("results/SL1344/leafcutter/fastqtl_output/naive_100kb_permuted.txt.gz") %>%
  dplyr::rename(transcript_id = gene_id)
ifng_fqtl = importFastQTLTable("results/SL1344/leafcutter/fastqtl_output/IFNg_100kb_permuted.txt.gz") %>%
  dplyr::rename(transcript_id = gene_id)
sl1344_fqtl = importFastQTLTable("results/SL1344/leafcutter/fastqtl_output/SL1344_100kb_permuted.txt.gz") %>%
  dplyr::rename(transcript_id = gene_id)
ifng_sl1344_fqtl = importFastQTLTable("results/SL1344/leafcutter/fastqtl_output/IFNg_SL1344_100kb_permuted.txt.gz") %>%
  dplyr::rename(transcript_id = gene_id)

fastqtl_pvalue_list = list(naive = naive_fqtl,
                           IFNg = ifng_fqtl,
                           SL1344 = sl1344_fqtl, 
                           IFNg_SL1344 = ifng_sl1344_fqtl)

#Remove duplicate transcript entries (do not know why they appeared)
fastqtl_pvalue_unique = purrr::map(fastqtl_pvalue_list, ~dplyr::group_by(.,transcript_id) %>% 
                                     dplyr::filter(row_number() == 1) %>% ungroup())

#Add transcript metadata to QTLs
fastqtl_pvalue_meta = purrr::map(fastqtl_pvalue_unique, ~dplyr::left_join(.,event_metadata_names, by = "transcript_id"))
saveRDS(fastqtl_pvalue_meta, "results/SL1344/leafcutter/leafcutter_min_pvalues_meta.rds")

#Apply bonferroni correction for p-values within cluster
fastqtl_bonferroni = purrr::map(fastqtl_pvalue_meta, ~dplyr::group_by(., cluster_id) %>% 
                                  dplyr::mutate(n_transcripts = length(transcript_id)) %>% 
                                  dplyr::arrange(cluster_id, p_beta) %>% dplyr::filter(row_number() == 1) %>% 
                                  dplyr::ungroup() %>% dplyr::mutate(p_bonferroni = p_beta * n_transcripts) %>% 
                                  dplyr::mutate(p_bonferroni = pmin(p_bonferroni, 1)) %>% 
                                  dplyr::mutate(p_fdr = p.adjust(p_bonferroni, method = "fdr")))
saveRDS(fastqtl_bonferroni,"results/SL1344/leafcutter/leafcutter_cluster_min_pvalues.rds")

#Call significant QTLs
leafcutter_qtl_hits = purrr::map_df(fastqtl_bonferroni, ~dplyr::filter(.,p_fdr < 0.1), .id = "condition_name")



