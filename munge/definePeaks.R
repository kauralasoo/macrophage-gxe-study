library("rtracklayer")
library("devtools")
library("dplyr")
library("rtracklayer")
load_all("../seqUtils/")

#### Work with H3K27Ac data ####
cond_A_names = c("H3K27Ac_rep1_A","H3K27Ac_rep2_A")
cond_A_narrowPeaks = loadNarrowPeaks("processed/Ivashkiv/", cond_A_names) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 
cond_A_broadPeaks = loadBroadPeaks("processed/Ivashkiv/", cond_A_names) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 

#Keep only the broad peaks that overlap a shared narrow peak
olaps = findOverlaps(cond_A_narrowPeaks, cond_A_broadPeaks)
cond_A_hc_peaks = cond_A_broadPeaks[unique(subjectHits(olaps)),]

cond_B_names = c("H3K27Ac_rep1_B","H3K27Ac_rep2_B")
cond_B_narrowPeaks = loadNarrowPeaks("processed/Ivashkiv/", cond_B_names) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 
cond_B_broadPeaks = loadBroadPeaks("processed/Ivashkiv/", cond_B_names) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 

#Keep only the broad peaks that overlap a shared narrow peak
olaps = findOverlaps(cond_B_narrowPeaks, cond_B_broadPeaks)
cond_B_hc_peaks = cond_B_broadPeaks[unique(subjectHits(olaps)),]

cond_C_names = c("H3K27Ac_rep1_C","H3K27Ac_rep2_C")
cond_C_narrowPeaks = loadNarrowPeaks("processed/Ivashkiv/", cond_C_names) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 
cond_C_broadPeaks = loadBroadPeaks("processed/Ivashkiv/", cond_C_names) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 

#Keep only the broad peaks that overlap a shared narrow peak
olaps = findOverlaps(cond_C_narrowPeaks, cond_C_broadPeaks)
cond_C_hc_peaks = cond_C_broadPeaks[unique(subjectHits(olaps)),]

cond_D_names = c("H3K27Ac_rep1_D","H3K27Ac_rep2_D")
cond_D_narrowPeaks = loadNarrowPeaks("processed/Ivashkiv/", cond_D_names) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 
cond_D_broadPeaks = loadBroadPeaks("processed/Ivashkiv/", cond_D_names) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 

#Keep only the broad peaks that overlap a shared narrow peak
olaps = findOverlaps(cond_D_narrowPeaks, cond_D_broadPeaks)
cond_D_hc_peaks = cond_D_broadPeaks[unique(subjectHits(olaps)),]

#Join peaks across all conditions together
all_peaks = listUnion(list(cond_A_hc_peaks, cond_B_hc_peaks, cond_C_hc_peaks, cond_D_hc_peaks))
metadata = data_frame(type = "exon", gene_id = paste("H3K27Ac_peak_", c(1:length(all_peaks)), sep = ""))
elementMetadata(all_peaks) = metadata

rtracklayer::export.bed(all_peaks, "annotations/H3K27Ac_joint_peaks.bed")
rtracklayer::export.gff3(all_peaks, "annotations/H3K27Ac_joint_peaks.gff3")

#Process STAT1 peaks
stat1_B_peaks = c("STAT1_rep1_B","STAT1_rep2_B")
stat1_B_narrowPeaks = loadNarrowPeaks("results/Ivashkiv/peak_calls", stat1_B_peaks, sub_dir = FALSE) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 

stat1_D_peaks = c("STAT1_rep1_D","STAT1_rep2_D")
stat1_D_narrowPeaks = loadNarrowPeaks("results/Ivashkiv/peak_calls", stat1_D_peaks, sub_dir = FALSE) %>% filterOverlaps(minOverlapCount = 2) %>% listUnion() 
stat1_union_peaks = union(stat1_D_narrowPeaks, stat1_B_narrowPeaks)
stat1_union_peaks = keepSeqlevels(stat1_union_peaks, c(1:22, "X","Y"))
stat1_metadata = data_frame(type = "exon", gene_id = paste("STAT1_peak_", c(1:length(stat1_union_peaks)), sep = ""))
elementMetadata(stat1_union_peaks) = stat1_metadata
  
rtracklayer::export.bed(stat1_union_peaks, "results/Ivashkiv/peak_calls/STAT1_joint_peaks.bed")
rtracklayer::export.gff3(stat1_union_peaks, "results/Ivashkiv/peak_calls/STAT1_joint_peaks.gff3")

