library("rtracklayer")
library("devtools")
library("dplyr")
library("rtracklayer")
load_all("../seqUtils/")

#CEBPbeta control
cebpb_ctrl_names = c("CEBPbeta_ctrl_201","CEBPbeta_ctrl_202","CEBPbeta_ctrl_203","CEBPbeta_ctrl_204")
CEBPbeta_ctrl_peaks = loadNarrowPeaks("results/OCallaghan/peak_calls/", cebpb_ctrl_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
rtracklayer::export.bed(CEBPbeta_ctrl_peaks, "results/OCallaghan/peak_calls/CEBPbeta_ctrl_consensus.bed")

#CEBPbeta oxLDL
cebpb_oxLDL_names = c("CEBPbeta_oxLDL_205","CEBPbeta_oxLDL_206","CEBPbeta_oxLDL_207","CEBPbeta_oxLDL_208")
CEBPbeta_oxLDL_peaks = loadNarrowPeaks("results/OCallaghan/peak_calls/", cebpb_oxLDL_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
rtracklayer::export.bed(CEBPbeta_oxLDL_peaks, "results/OCallaghan/peak_calls/CEBPbeta_oxLDL_consensus.bed")

#H3K27Ac ctrl
H3K27Ac_ctrl_names = c("H3K27ac_ctrl_250", "H3K27ac_ctrl_252")
H3K27Ac_ctrl_peaks = loadNarrowPeaks("results/OCallaghan/peak_calls/", H3K27Ac_ctrl_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
H3K27Ac_ctrl_broadPeaks = loadBroadPeaks("results/OCallaghan/peak_calls/", H3K27Ac_ctrl_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()

#Keep only the broad peaks that overlap a shared narrow peak
olaps = findOverlaps(H3K27Ac_ctrl_peaks, H3K27Ac_ctrl_broadPeaks)
H3K27Ac_ctrl_hc_peaks = H3K27Ac_ctrl_broadPeaks[unique(subjectHits(olaps)),]
rtracklayer::export.bed(H3K27Ac_ctrl_hc_peaks, "results/OCallaghan/peak_calls/H2K27Ac_ctrl_consensus.bed")

#H3K27Ac oxLDL
H3K27Ac_oxLDL_names = c("H3K27ac_oxLDL_249", "H3K27ac_oxLDL_251")
H3K27Ac_oxLDL_peaks = loadNarrowPeaks("results/OCallaghan/peak_calls/", H3K27Ac_oxLDL_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
H3K27Ac_oxLDL_broadPeaks = loadBroadPeaks("results/OCallaghan/peak_calls/", H3K27Ac_oxLDL_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()

#Keep only the broad peaks that overlap a shared narrow peak
olaps = findOverlaps(H3K27Ac_oxLDL_peaks, H3K27Ac_oxLDL_broadPeaks)
H3K27Ac_oxLDL_hc_peaks = H3K27Ac_oxLDL_broadPeaks[unique(subjectHits(olaps)),]
rtracklayer::export.bed(H3K27Ac_oxLDL_hc_peaks, "results/OCallaghan/peak_calls/H2K27Ac_oxLDL_consensus.bed")

