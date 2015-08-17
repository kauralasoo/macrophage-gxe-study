library("rtracklayer")
library("devtools")
library("dplyr")
load_all("../macrophage-gxe-study/macrophage-gxe-study/seqUtils/")

#Import sample names
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

a = loadCounts("processed/SL1344/", sample_names[1:5], counts_suffix = ".consensus_peaks.counts.txt")

heatmap.2(cor(a[,3:ncol(a)],method = "spearman"))