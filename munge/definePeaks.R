library("rtracklayer")
library("devtools")
library("dplyr")
load_all("../seqUtils/")
load_all("macrophage-chromatin/housekeeping//")

#Load all sample names from disk
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
design = housekeeping::constructDesignMatrix_ATAC(sample_names)

#Find peaks in condition A
cond_A_names = dplyr::filter(design, condition_char == "A")$sample_id
cond_A_peaks = loadNarrowPeaks("processed/SL1344/", cond_A_names)
cond_A_peak_count = lapply(cond_A_peaks, length) %>% 
  plyr::ldply(.id = "sample_id") %>% 
  dplyr::rename(peak_count = V1) 
cond_A_union_all = listUnion(cond_A_peaks)

#Filter peaks by overlap
cond_A_filtered = filterOverlaps(cond_A_peaks, 2)
cond_A_union = listUnion(cond_A_filtered)

#Find peaks in condition B
cond_B_names = dplyr::filter(design, condition_char == "B")$sample_id
cond_B_peaks = loadNarrowPeaks("processed/SL1344/", cond_B_names)
cond_B_peak_count = lapply(cond_B_peaks, length) %>% 
  plyr::ldply(.id = "sample_id") %>% 
  dplyr::rename(peak_count = V1)
cond_B_filtered = filterOverlaps(cond_B_peaks, 2)
cond_B_union = listUnion(cond_B_filtered)

#Find peaks in condition C
cond_C_names = dplyr::filter(design, condition_char == "C")$sample_id
cond_C_peaks = loadNarrowPeaks("processed/SL1344/", cond_C_names)
cond_C_peak_count = lapply(cond_C_peaks, length) %>% 
  plyr::ldply(.id = "sample_id") %>% 
  dplyr::rename(peak_count = V1)
cond_C_filtered = filterOverlaps(cond_C_peaks, 2)
cond_C_union = listUnion(cond_C_filtered)

#Find peaks in condition D
cond_D_names = dplyr::filter(design, condition_char == "D")$sample_id
cond_D_peaks = loadNarrowPeaks("processed/SL1344/", cond_D_names)
cond_D_peak_count = lapply(cond_D_peaks, length) %>% 
  plyr::ldply(.id = "sample_id") %>% 
  dplyr::rename(peak_count = V1)
cond_D_filtered = filterOverlaps(cond_D_peaks, 2)
cond_D_union = listUnion(cond_D_filtered)

#Join all peaks together
all_peaks = listUnion(list(cond_A_union, cond_B_union, cond_C_union, cond_D_union))
metadata = data_frame(type = "exon", gene_id = paste("ATAC_peak_", c(1:length(all_peaks)), sep = ""))
elementMetadata(all_peaks) = metadata

#Export as a gff file
rtracklayer::export.gff3(all_peaks, "annotations/ATAC_Seq_joint_peaks.gff3")
rtracklayer::export.bed(all_peaks, "annotations/ATAC_Seq_joint_peaks.bed")

#Use peak count per sample as QC measures
peak_counts = rbind(cond_A_peak_count, cond_B_peak_count, cond_C_peak_count, cond_D_peak_count)
write.table(peak_counts, "macrophage-chromatin/data/SL1344/QC_measures/macs2_peaks_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)
