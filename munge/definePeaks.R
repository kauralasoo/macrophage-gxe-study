library("rtracklayer")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Find peaks in condition A
cond_A_names = c("bima_A_ATAC","vass_A_ATAC","cicb_A_ATAC","eofe_A_ATAC")
cond_A_peaks = loadNarrowPeaks("processed/SL1344/", cond_A_names)
cond_A_filtered = filterOverlaps(cond_A_peaks, 3)
cond_A_union = listUnion(cond_A_filtered)

#Find peaks in condition B
cond_B_names = c("bima_B_ATAC","vass_B_ATAC","cicb_B_ATAC","eofe_B_ATAC")
cond_B_peaks = loadNarrowPeaks("processed/SL1344/", cond_B_names)
cond_B_filtered = filterOverlaps(cond_B_peaks, 3)
cond_B_union = listUnion(cond_B_filtered)

#Find peaks in condition C
cond_C_names = c("bima_C_ATAC","vass_C_ATAC","cicb_C_ATAC","eofe_C_ATAC")
cond_C_peaks = loadNarrowPeaks("processed/SL1344/", cond_C_names)
cond_C_filtered = filterOverlaps(cond_C_peaks, 3)
cond_C_union = listUnion(cond_C_filtered)

#Find peaks in condition D
cond_D_names = c("bima_D_ATAC","vass_D_ATAC","cicb_D_ATAC","eofe_D_ATAC")
cond_D_peaks = loadNarrowPeaks("processed/SL1344/", cond_D_names)
cond_D_filtered = filterOverlaps(cond_D_peaks, 3)
cond_D_union = listUnion(cond_D_filtered)

#Join all peaks together
all_peaks = listUnion(list(cond_A_union, cond_B_union, cond_C_union, cond_D_union))
metadata = data_frame(type = "exon", gene_id = paste("ATAC_peak_", c(1:length(all_peaks)), sep = ""))
elementMetadata(all_peaks) = metadata

#Export as a gff file
rtracklayer::export.gff3(all_peaks, "annotations/ATAC_Seq_joint_peaks.gff3")
rtracklayer::export.bed(all_peaks, "annotations/ATAC_Seq_joint_peaks.bed")
