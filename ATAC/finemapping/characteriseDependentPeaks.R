library("dplyr")
library("GenomicRanges")
library("devtools")
library("rtracklayer")
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
master_peak_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% all_dependent_peaks$master_id) %>% 
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

#Find dependent overlaps
dep_olaps = findOverlaps(dependent_peak_ranges, tx_ranges)
dep_peaks = elementMetadata(dependent_peak_ranges[queryHits(dep_olaps)])
dep_transcripts = elementMetadata(tx_ranges[subjectHits(dep_olaps)])
dep_matches = cbind(dep_peaks, dep_transcripts) %>% tbl_df2()

#Filter overlaps
dep_filtered_overlaps = dplyr::filter(dep_matches, transcript_biotype %in% c("protein_coding", "lincRNA")) %>% 
  dplyr::select(peak_id, gene_id, gene_name) %>% unique()

#Find master overlaps
master_olaps = findOverlaps(master_peak_ranges, tx_ranges)
master_peaks = elementMetadata(master_peak_ranges[queryHits(master_olaps)])
master_transcripts = elementMetadata(tx_ranges[subjectHits(master_olaps)])
master_matches = cbind(master_peaks, master_transcripts) %>% tbl_df2()

#Filter overlaps
master_filtered_overlaps = dplyr::filter(master_matches, transcript_biotype %in% c("protein_coding", "lincRNA")) %>% 
  dplyr::select(peak_id, gene_id, gene_name) %>% unique()

##### Count overlaps with h2k27ac peaks ####
h3k27ac_peaks = import.gff3("annotations/chromatin/H3K27Ac_joint_peaks.gff3")

olaps = findOverlaps(master_peak_ranges, h3k27ac_peaks)
master_peaks = elementMetadata(master_peak_ranges[queryHits(olaps)]) %>% tbl_df2() %>% unique()

olaps = findOverlaps(dependent_peak_ranges, h3k27ac_peaks)
dependent_peaks = elementMetadata(dependent_peak_ranges[queryHits(olaps)]) %>% tbl_df2() %>% unique()


#Combine results into a table
result_df = data_frame(type = c("master", "dependent"), 
                       total_count = c(length(master_peak_ranges), length(dependent_peak_ranges)),
                       promoter_count = c(nrow(master_filtered_overlaps), nrow(dep_filtered_overlaps)),
                       H3K27ac_count = c(nrow(master_peaks), nrow(dependent_peaks))) %>%
  dplyr::mutate(promoter_fraction = round(promoter_count/total_count,3),
                H3K27ac_fraction = round(H3K27ac_count/total_count,3))
write.table(result_df, "figures/tables/enhancer_promoter_overlap.txt", sep = "\t", quote = FALSE, row.names = FALSE)







