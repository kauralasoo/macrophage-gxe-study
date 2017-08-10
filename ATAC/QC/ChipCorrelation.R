library("devtools")
library("dplyr")
library("ggplot2")
library("purrr")
library("GenomicFeatures")
load_all("../seqUtils/")
load_all("../wiggleplotr")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#ATAC-seq
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, 
                                                  levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
atac_granges = dplyr::transmute(atac_list$gene_metadata, seqnames = chr, start, end, strand, gene_id) %>% dataFrameToGRanges()

#Import PU1 counts
sample_names = read.table("macrophage-gxe-study/data/chromatin/ChIP/Schultze_sample_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
schultze_counts = loadCounts("processed/Schultze/counts/", sample_names, counts_suffix = ".ATAC_peaks.counts.txt", sub_dir = FALSE)
schultze_matrix = schultze_counts[,sample_names]
row.names(schultze_matrix) = schultze_counts$gene_id

#Load Schultze dataset
schultze_names = c("PU1_naive")
schultze_peaks = loadNarrowPeaks("processed/Schultze/peak_calls/", schultze_names, sub_dir = FALSE)

overlap_atac = atac_granges[queryHits(findOverlaps(atac_granges, schultze_peaks$PU1_naive))]
overlap_peaks = overlap_atac$gene_id

atac = atac_list$counts[overlap_peaks,]
chip = schultze_matrix[overlap_peaks,]
pdf("figures/supplementary/ATAC_vs_PU1_scatter.pdf", width = 5, height = 5)
smoothScatter(log(atac$bima_A_ATAC + 1, 2), log(chip$PU1_naive +1, 2), xlab = "Log2(ATAC_signal + 1)", ylab = "Log2(PU1_signal + 1)")
dev.off()
cor(log(atac$bima_A_ATAC + 1, 2), log(chip$PU1_naive +1, 2))


#Import CEBPb counts
sample_names = read.table("macrophage-gxe-study/data/chromatin/ChIP/OCallaghan_sample_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
ocallaghan_counts = loadCounts("processed/OCallaghan/counts/", sample_names, counts_suffix = ".ATAC_peaks.counts.txt", sub_dir = FALSE)
ocallaghan_matrix = ocallaghan_counts[,sample_names]
row.names(ocallaghan_matrix) = ocallaghan_counts$gene_id

ocallaghan_peaks = loadNarrowPeaks("processed/OCallaghan/peak_calls/", sample_names, sub_dir = FALSE)

overlap_atac = atac_granges[queryHits(findOverlaps(atac_granges, ocallaghan_peaks$CEBPbeta_ctrl_201))]
overlap_peaks = overlap_atac$gene_id

atac = atac_data$counts[overlap_peaks,]
chip = ocallaghan_matrix[overlap_peaks,]
smoothScatter(log(atac$bima_A_ATAC + 1, 2), log(chip$CEBPbeta_ctrl_201 +1, 2))
cor(log(atac$bima_A_ATAC + 1, 2), log(chip$CEBPbeta_ctrl_201 +1, 2))







