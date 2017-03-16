library("dplyr")
library("GenomicRanges")
library("devtools")
load_all("../seqUtils/")

#Import transcript annotations
tx_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_79/Homo_sapiens.GRCh38.79.transcript_data.rds") %>%
  dplyr::tbl_df()

#Import peak types
result_list = readRDS("results/ATAC/QTLs/qtl_peak_type_assignment.rds")
all_dependent_peaks = purrr::map_df(result_list$dependents, ~dplyr::select(.,master_id, dependent_id)) %>% 
  unique()

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
dependent_peak_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% all_dependent_peaks$dependent_id) %>% 
  dplyr::transmute(seqnames = chr, start, end, strand, peak_id = gene_id) %>% 
  dataFrameToGRanges()

tx_ranges = dplyr::filter(tx_metadata, gene_biotype %in% c("protein_coding", "lincRNA")) %>% 
  dplyr::transmute(seqnames = chromosome_name, 
                   start = ifelse(strand == 1, transcript_start, transcript_end), 
                   gene_id = ensembl_gene_id, 
                   transcript_id = ensembl_transcript_id,
                   strand = "+", transcript_biotype,
                   gene_name = external_gene_name) %>% 
  dplyr::mutate(end = start) %>%
  dplyr::mutate(start = start - 500, end = end + 500) %>%
  dataFrameToGRanges()

olaps = findOverlaps(dependent_peak_ranges, tx_ranges)
peaks = elementMetadata(dependent_peak_ranges[queryHits(olaps)])
transcripts = elementMetadata(tx_ranges[subjectHits(olaps)])
matches = cbind(peaks, transcripts) %>% tbl_df2()

#Filter overlaps
filtered_overlaps = dplyr::filter(matches, transcript_biotype %in% c("protein_coding", "lincRNA")) %>% 
  dplyr::select(peak_id, gene_id, gene_name) %>% unique()

#Import condition_specific eQTLs
variable_qtls = readRDS("results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")
gene_clusters = dplyr::select(variable_qtls$appear, gene_id, snp_id, max_condition) %>% ungroup() %>% unique()

dplyr::semi_join(filtered_overlaps, gene_clusters)
