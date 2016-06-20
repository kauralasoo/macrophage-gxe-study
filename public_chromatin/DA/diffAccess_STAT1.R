library("dplyr")
library("DESeq2")
library("devtools")
library("edgeR")
library("rtracklayer")
load_all("macrophage-chromatin/housekeeping/")
load_all("../seqUtils/")

#Import count data
stat1_counts = readRDS("results/Ivashkiv/STAT1_combined_counts.rds")
sample_names =  c("STAT1_rep1_B","STAT1_rep1_D","STAT1_rep2_B","STAT1_rep2_D")
stat1_design = constructDesignMatrix_Ivashkiv(sample_names)

length_df = dplyr::select(stat1_counts, gene_id, length)
rownames(stat1_counts) = stat1_counts$gene_id
diff_counts = stat1_counts[,stat1_design$sample_id]

dge = DGEList(counts = diff_counts)
dge = calcNormFactors(dge, method = "TMM")

#Apply voom transformation
design_matrix = model.matrix(~condition_name, stat1_design)
v <- voom(dge,design_matrix,plot=TRUE)
fit = lmFit(v, design_matrix)
fit <- eBayes(fit)

#Identify regulated peaks
stat1_hits = topTable(fit, coef = c(2),  lfc = 1, p.value = 0.1, number = 1e5) %>% tidyTopTable()
stat1_syn_up = dplyr::filter(stat1_hits, logFC > 0)
stat1_syn_down = dplyr::filter(stat1_hits, logFC < 0)

#Make a bed file of STAT1 peaks
stat1_peaks = import.gff3("results/Ivashkiv/peak_calls/STAT1_joint_peaks.gff3")
stat1_syn_up_peaks = stat1_peaks[stat1_peaks$gene_id %in% stat1_syn_up$gene_id,]
stat1_syn_down_peaks = stat1_peaks[stat1_peaks$gene_id %in% stat1_syn_down$gene_id,]

#Add names
elementMetadata(stat1_peaks)$name = "STAT1"
elementMetadata(stat1_syn_up_peaks)$name = "STAT1_syn_up"
elementMetadata(stat1_syn_down_peaks)$name = "STAT1_syn_down"

#Export BED
stat1_all_peaks = c(stat1_peaks, stat1_syn_up_peaks, stat1_syn_down_peaks)
elementMetadata(stat1_all_peaks) = elementMetadata(stat1_all_peaks)[,-c(3,4)]
export.bed(stat1_all_peaks, "results/Ivashkiv/DA/STAT1_grouped_peaks.bed")



