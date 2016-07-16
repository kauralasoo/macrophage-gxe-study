library("rtracklayer")
library("devtools")
library("dplyr")
library("rtracklayer")
load_all("../seqUtils/")

#Import sample names
sample_names = as.vector(read.table("macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt")[,1])
sample_names = sample_names[!(sample_names %in% c("BTS0010_input","BTS0011_input"))]

#Import peak calls
irf2_narrowPeaks = loadNarrowPeaks("/Volumes/JetDrive/databases/peak_calls/Fairfax_et_al/", 
                                   sample_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion() 
saveRDS(irf2_narrowPeaks, "results/public_chromatin/joint_peaks/IRF2_Fairfax_peaks.rds")


