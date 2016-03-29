library("rtracklayer")
library("devtools")
library("dplyr")
library("rtracklayer")
load_all("../seqUtils/")

#naive CIITA
ciita_naive_names = c("MO_rep1_naive_CIITA","MO_rep2_naive_CIITA")
ciita_naive_peaks = loadNarrowPeaks("results/Knight/peak_calls/", ciita_naive_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
elementMetadata(ciita_naive_peaks)$name = "CIITA_naive"

#IFNg CIITA
ciita_ifng_names = c("MO_rep1_IFNg_CIITA","MO_rep2_IFNg_CIITA")
ciita_ifng_peaks = loadNarrowPeaks("results/Knight/peak_calls/", ciita_ifng_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
elementMetadata(ciita_ifng_peaks)$name = "CIITA_IFNg"

#naive RFX5
rfx5_naive_names = c("MO_rep1_naive_RFX5","MO_rep2_naive_RFX5")
rfx5_naive_peaks = loadNarrowPeaks("results/Knight/peak_calls/", rfx5_naive_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
elementMetadata(rfx5_naive_peaks)$name = "RFX5_naive"


rfx5_ifng_names = c("MO_rep1_IFNg_RFX5","MO_rep2_IFNg_RFX5")
rfx5_ifng_peaks = loadNarrowPeaks("results/Knight/peak_calls/", rfx5_ifng_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
elementMetadata(rfx5_ifng_peaks)$name = "RFX5_IFNg"

joint_peaks = c(ciita_naive_peaks, ciita_ifng_peaks, rfx5_naive_peaks, rfx5_ifng_peaks)
export.bed(joint_peaks, "results/Knight/CIITA-RFX5_joint_peaks.bed")

