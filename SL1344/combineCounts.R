library("devtools")
library("dplyr")
load_all("macrophage-gxe-study/seqUtils/")

sample_names = read.table("fastq/SL1344_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

#Import Gencode Basic counts
basic_data = loadCounts("STAR/SL1344/", sample_names, sub_dir = TRUE, counts_suffix = ".gencode_basic.counts.txt")
write.table(basic_data, "results/SL1344/SL1344_basic_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Import complete Ensembl 79 counts
complete_data = loadCounts("STAR/SL1344/", sample_names, sub_dir = TRUE, counts_suffix = ".counts.txt")
write.table(complete_data, "results/SL1344/SL1344_complete_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

