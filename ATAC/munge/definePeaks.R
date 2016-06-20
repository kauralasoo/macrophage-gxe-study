library("rtracklayer")
library("devtools")
library("dplyr")
load_all("../seqUtils/")
load_all("macrophage-chromatin/housekeeping//")

#Load all sample names from disk
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
design = housekeeping::constructDesignMatrix_ATAC(sample_names)

#Import peak call from each condition
#Find peaks in condition A
cond_A_names = dplyr::filter(design, condition_char == "A")$sample_id
cond_A_peaks = loadNarrowPeaks("processed/SL1344/", cond_A_names)
cond_A_peak_count = lapply(cond_A_peaks, length) %>% 
  plyr::ldply(.id = "sample_id") %>% 
  dplyr::rename(peak_count = V1) 
cond_A_union_all = listUnion(cond_A_peaks)

#Find peaks in condition B
cond_B_names = dplyr::filter(design, condition_char == "B")$sample_id
cond_B_peaks = loadNarrowPeaks("processed/SL1344/", cond_B_names)
cond_B_peak_count = lapply(cond_B_peaks, length) %>% 
  plyr::ldply(.id = "sample_id") %>% 
  dplyr::rename(peak_count = V1)

#Find peaks in condition C
cond_C_names = dplyr::filter(design, condition_char == "C")$sample_id
cond_C_peaks = loadNarrowPeaks("processed/SL1344/", cond_C_names)
cond_C_peak_count = lapply(cond_C_peaks, length) %>% 
  plyr::ldply(.id = "sample_id") %>% 
  dplyr::rename(peak_count = V1)

#Find peaks in condition D
cond_D_names = dplyr::filter(design, condition_char == "D")$sample_id
cond_D_peaks = loadNarrowPeaks("processed/SL1344/", cond_D_names)
cond_D_peak_count = lapply(cond_D_peaks, length) %>% 
  plyr::ldply(.id = "sample_id") %>% 
  dplyr::rename(peak_count = V1)

#Make peak list
peak_list = list(naive = cond_A_peaks, IFNg = cond_B_peaks, SL1344 = cond_C_peaks, IFNg_SL1344 = cond_D_peaks)
saveRDS(peak_list, "results/ATAC/ATAC_peak_list.rds")
peak_list = readRDS("results/ATAC/ATAC_peak_list.rds")



#Filter peaks by overlap
naive_filtered = filterOverlaps(peak_list$naive, 3)
naive_union = listUnion(naive_filtered)

ifng_peaks = peak_list$IFNg[!(names(peak_list$IFNg) == "fikt_B_ATAC")]
ifng_peaks_filtered = filterOverlaps(ifng_peaks, 3)
IFNg_union = listUnion(ifng_peaks_filtered)

SL1344_filtered = filterOverlaps(peak_list$SL1344, 3)
SL1344_union = listUnion(SL1344_filtered)

IFNg_SL1344_filtered = filterOverlaps(peak_list$IFNg_SL1344, 3)
IFNg_SL1344_union = listUnion(IFNg_SL1344_filtered)


#Join all peaks together
all_peaks = listUnion(list(naive_union, IFNg_union, SL1344_union, IFNg_SL1344_union))
metadata = data_frame(type = "exon", gene_id = paste("ATAC_peak_", c(1:length(all_peaks)), sep = ""))
elementMetadata(all_peaks) = metadata

#Add strand
strand(all_peaks) = "+"

#Export as a gff file
rtracklayer::export.gff3(all_peaks, "annotations/ATAC_consensus_peaks.gff3")
rtracklayer::export.bed(all_peaks, "annotations/ATAC_consensus_peaks.bed")

#Use peak count per sample as QC measures
peak_counts = rbind(cond_A_peak_count, cond_B_peak_count, cond_C_peak_count, cond_D_peak_count)
write.table(peak_counts, "macrophage-chromatin/data/SL1344/QC_measures/macs2_peaks_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)
