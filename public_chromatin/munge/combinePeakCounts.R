library("rtracklayer")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import sample names
sample_names = read.table("macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
sample_names = sample_names[grepl("H3K27Ac",sample_names)]

h3k27ac_counts = loadCounts("processed/Ivashkiv//", sample_names, counts_suffix = ".H3K27Ac_joint_peaks.counts.txt")
write.table(h3k27ac_counts, "results/Ivashkiv/H3K27Ac_combined_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(h3k27ac_counts, "results/Ivashkiv/H3K27Ac_combined_counts.rds")

#IRF1
irf_names = c("IRF1_A","IRF1_B","IRF1_E","IRF1_F")
irf_counts = loadCounts("processed/Ivashkiv/", irf_names, counts_suffix = ".IRF1_joint_peaks.counts.txt")
write.table(irf_counts, "results/Ivashkiv/IRF1_combined_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(irf_counts, "results/Ivashkiv/IRF1_combined_counts.rds")

#STAT1
stat1_names = c("STAT1_rep1_A","STAT1_rep1_B","STAT1_rep1_C","STAT1_rep1_D","STAT1_rep2_B","STAT1_rep2_D")
stat_counts = loadCounts("processed/Ivashkiv/", stat1_names, counts_suffix = ".STAT1_joint_peaks.counts.txt")
write.table(stat_counts, "results/Ivashkiv/STAT1_combined_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(stat_counts, "results/Ivashkiv/STAT1_combined_counts.rds")
