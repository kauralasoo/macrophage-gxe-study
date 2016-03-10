library("devtools")
load_all("../seqUtils/")

processed/Rehli/MAC_PU1/MAC_PU1_peaks.narrowPeak

#Load Rehli dataset
sample_names = c("MAC_PU1","MAC_CEBPbeta","CTCF_MAC")
rehli_peaks = loadNarrowPeaks("processed/Rehli/", sample_names)

#Load Ivashkiv peaks
sample_names = c("STAT1_rep1_A","STAT1_rep1_B","STAT1_rep1_C","STAT1_rep1_D","STAT1_rep2_B","STAT1_rep2_D",
                 "IRF1_A","IRF1_B","IRF1_E","IRF1_F")
ivashkiv_peaks = loadNarrowPeaks("processed/Ivashkiv/", sample_names)

peak_list = list(Rehli = rehli_peaks, Ivashkiv = ivashkiv_peaks)
saveRDS(peak_list, "results/ATAC/DA/Chip_peak_lists.rds")

