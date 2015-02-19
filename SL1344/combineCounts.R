library("DESeq2")
library("dplyr")
load_all("macrophage-gxe-study/seqUtils/")

sample_names = read.table("fastq/SL1344_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

data = loadCounts("STAR/SL1344/", sample_names)
write.table(data, "results/SL1344/combined_counts.txt", sep = "\t", quote = FALSE)
