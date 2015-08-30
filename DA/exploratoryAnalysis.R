library(gplots)
source("macrophage-chromatin/housekeeping/constructDesignMatrices.r")

#Import peak counts
atac_counts = readRDS("results/ATAC/ATAC_combined_counts.rds")

#Construct design matrix
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
design = constructDesignMatrix_ATAC(sample_names)

#Make a counts matrix
counts_df = atac_counts[,-c(1,2)]
rownames(counts_df) = atac_counts$gene_id
lengths_df = atac_counts[,1:2]

#Calculate TPM values
atac_tpm = calculateTPM(counts_df, lengths_df)
atac_tpm_log2 = log(atac_tpm + 0.1, 2)

#Make heatmap of TPM counts
pdf("results/ATAC/QC/ATAC_spearman_heatmap.pdf", width = 10, height = 10)
heatmap.2(cor(atac_tpm_log2,method = "spearman"), margins = c(10,10), tracecol = NA)
dev.off()

#Make a PCA plot of log2 TPM values
pca = performPCA(atac_tpm_log2, design)
pca_plot = ggplot(pca$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + 
  geom_text()
ggsave("results/ATAC/QC/ATAC_tpm_pca.pdf", pca_plot, width = 7, height = 6)
