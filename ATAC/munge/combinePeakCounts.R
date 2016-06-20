library("rtracklayer")
library("devtools")
library("dplyr")
library("gplots")
load_all("../seqUtils/")

#Import sample names
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
atac_counts = loadCounts("processed/SL1344/", sample_names, counts_suffix = ".consensus_peaks.counts.txt")

#Make heatmap
pdf("results/ATAC/QC/ATAC_spearman_heatmap.pdf", width = 10, height = 10)
heatmap.2(cor(atac_counts[,3:ncol(atac_counts)],method = "spearman"), margins = c(10,10), tracecol = NA)
dev.off()

#Save counts to disk
write.table(atac_counts, "results/ATAC/ATAC_combined_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(atac_counts, "results/ATAC/ATAC_combined_counts.rds")
