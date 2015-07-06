library("DESeq2")
library("dplyr")
library("devtools")
load_all("macrophage-gxe-study/seqUtils/")

sample_names = read.table("macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_all.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

#Import GENCODE Basic counts
data = loadCounts("STAR/acLDL/", sample_names, counts_suffix = ".gencode_basic.counts.txt")
write.table(data, "results/acLDL/acLDL_basic_counts.txt", sep = "\t", quote = FALSE)

#Import all counts
data_full = loadCounts("STAR/acLDL/", sample_names, counts_suffix = ".counts.txt")
write.table(data_full, "results/acLDL/acLDL_complete_counts.txt", sep = "\t", quote = FALSE)
