library(gplots)
library(ggplot2)
source("macrophage-chromatin/housekeeping/constructDesignMatrices.r")

#Import peak counts
h3k27ac_counts = readRDS("results/Ivashkiv/H3K27Ac_combined_counts.rds")

#Construct design matrix
sample_names = read.table("macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
sample_names = sample_names[grepl("H3K27Ac",sample_names)]
design = constructDesignMatrix_Ivashkiv(sample_names)

#Make a counts matrix
counts_df = h3k27ac_counts[,-c(1,2)]
rownames(counts_df) = h3k27ac_counts$gene_id
lengths_df = h3k27ac_counts[,1:2]

#Calculate TPM values
h3k27ac_tpm = calculateTPM(counts_df, lengths_df)
h3k27ac_tpm_log2 = log(h3k27ac_tpm + 0.1, 2)

#Make heatmap of TPM counts
pdf("results/Ivashkiv//QC/H3K27Ac_spearman_heatmap.pdf", width = 10, height = 10)
heatmap.2(cor(h3k27ac_tpm_log2,method = "pearson"), margins = c(10,10), tracecol = NA)
dev.off()

#Make a PCA plot of log2 TPM values
pca = performPCA(h3k27ac_tpm_log2, design)
pca_plot = ggplot(pca$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + 
  geom_text()
ggsave("results/Ivashkiv/QC/H3K27Ac_tpm_pca.pdf", pca_plot, width = 9, height = 8)