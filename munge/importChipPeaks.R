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


#### Make consensus peaks ####
chr_name = read.table("../macrophage-gxe-study/macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt", stringsAsFactors = FALSE)

#IRF1
irf_names = c("IRF1_A","IRF1_B","IRF1_E","IRF1_F")
irf_peaks = loadNarrowPeaks("processed/Ivashkiv/", irf_names)
irf_peaks = makeUnionPeaks(irf_peaks, chr_name[,1], "IRF1_peak_")
rtracklayer::export.gff3(irf_peaks, "annotations/IRF1_joint_peaks.gff3")

#STAT1
stat1_names = c("STAT1_rep1_A","STAT1_rep1_B","STAT1_rep1_C","STAT1_rep1_D")
stat1_peaks = loadNarrowPeaks("processed/Ivashkiv/", stat1_names)
stat1_peaks = makeUnionPeaks(stat1_peaks, chr_name[,1], "STAT1_peak_")
rtracklayer::export.gff3(stat1_peaks, "annotations/STAT1_joint_peaks.gff3")
