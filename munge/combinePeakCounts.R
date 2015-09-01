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