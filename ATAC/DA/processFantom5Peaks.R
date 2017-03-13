library("dplyr")
library("tidyr")
library("devtools")
library(GenomicRanges)
load_all("../seqUtils/")

#Import fantom ata
fantom_peaks = readr::read_tsv("databases/FANTOM5/hg19.cage_peak_phase1and2combined_counts_ann_decoded.osc.txt.gz.extract.tsv")
colnames(fantom_peaks) = c("peak_id", "Ctrl_1", "Ctrl_2", "Ctrl_3", "LPS_1", "LPS_2", "LPS_3")
fantom_peaks[,2:7] = t(t(fantom_peaks[,2:7])/colSums(fantom_peaks[,2:7]))*1e6

#Import atac data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Import caQTL and eQTL pairs
caQTL_eQTL_pairs = readRDS("results/ATAC_RNA_overlaps/caQTL_eQTL_pairs_betas.rds")
filtered_pairs = purrr::map_df(caQTL_eQTL_pairs, identity, .id = "condition") %>% 
  dplyr::filter(condition_name == "naive", phenotype == "ATAC-seq")

#Identify two groups of peaks
appear_peaks = dplyr::filter(filtered_pairs, condition == "SL1344", abs(beta) < 0.59) %>% 
  dplyr::select(peak_id) %>% unique()
persistent_peaks = dplyr::filter(filtered_pairs, condition == "SL1344", abs(beta) > 0.59) %>% 
  dplyr::select(peak_id) %>% unique()


#Extract coordinates and strand
cage_bed = tidyr::separate(fantom_peaks, peak_id, c("chr", "other"), sep = ":", remove = FALSE) %>%
  tidyr::separate(other, c("start","end"), sep = "\\.\\.") %>%
  tidyr::separate(end, c("end","strand"), sep = ",") %>%
  dplyr::mutate(score = 1) %>%
  dplyr::select(chr,start,end, peak_id, score, strand, everything())
write.table(cage_bed, "databases/FANTOM5/CAGE_peak_counts.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Import new peak coordinates (after CrossMap to GRCh38)
new_coords = readr::read_tsv("databases/FANTOM5/CAGE_peak_counts.hg38.bed", col_names = c("chr","start","end","peak_id")) %>%
  tidyr::separate(chr, c("none","chr"), sep = "chr") %>%
  dplyr::select(-none)

#Identify cage peaks that overlap ATAC peaks
peak_ranges = dplyr::transmute(new_coords, seqnames = chr, start, end, strand = "+", peak_id) %>%
  dataFrameToGRanges()
atac_ranges = dplyr::transmute(atac_list$gene_metadata, seqnames = chr, start, end, strand, gene_id) %>% 
  dataFrameToGRanges()
olaps = findOverlaps(atac_ranges, peak_ranges)
map = data_frame(peak_id = atac_ranges[queryHits(olaps)]$gene_id, 
                 fantom_peak_id = peak_ranges[subjectHits(olaps)]$peak_id)

#Condition-specifc overlaps
overlapping_peaks = dplyr::filter(map, peak_id %in% appear_peaks$peak_id)

a = dplyr::filter(fantom_peaks, peak_id %in% overlapping_peaks$fantom_peak_id) %>% 
  dplyr::rename(fantom_peak_id = peak_id) %>% dplyr::left_join(map) %>%
  dplyr::select(peak_id, Ctrl_1:LPS_3) %>%
  tidyr::gather("sample_id", "expression", Ctrl_1:LPS_3) %>% 
  group_by(peak_id, sample_id) %>%
  dplyr::summarise(expression = sum(expression)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(condition_name = ifelse(grepl("Ctrl", sample_id), "naive", "LPS")) %>%
  dplyr::mutate(log_exp = log(expression + 1,2)) %>%
  dplyr::group_by(peak_id, condition_name) %>%
  dplyr::mutate(mean_log_exp = mean(log_exp)) %>%
  dplyr::group_by(peak_id) %>% 
  dplyr::mutate(max_log_exp = max(mean_log_exp)) %>%
  dplyr::filter(max_log_exp > 1) %>%
  dplyr::ungroup()
ggplot(a, aes(x = condition_name, y = log_exp)) + geom_point() + facet_wrap(~peak_id)




