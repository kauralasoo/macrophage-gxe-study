library("dplyr")
library("pheatmap")
library("devtools")
library("ggplot2")
load_all("macrophage-chromatin/housekeeping/")
load_all("../seqUtils/")

#Import peak counts
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_data$design$condition_name = factor(atac_data$design$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))

#Construct design matrix
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
design = constructDesignMatrix_ATAC(sample_names)

#Make a counts matrix
counts_df = atac_counts[,-c(1,2)]
rownames(counts_df) = atac_counts$gene_id
lengths_df = atac_counts[,1:2]

#Calculate TPM values
atac_tpm = calculateTPM(atac_data$exprs_counts, atac_data$gene_metadata)
atac_tpm_log2 = log(atac_tpm + 0.1, 2)

#Make heatmap of TPM counts
pdf("results/ATAC/QC/ATAC_spearman_heatmap.pdf", width = 10, height = 10)
heatmap.2(cor(atac_tpm_log2,method = "spearman"), margins = c(10,10), tracecol = NA)
dev.off()

#Make a heatmap of CQN normalized accesibility data
meta = as.data.frame(dplyr::select(atac_data$design, condition_name))
rownames(meta) = atac_data$design$sample_id
pheatmap(cor(atac_tpm_log2, method = "spearman"), annotation_row = meta, 
         show_rownames = FALSE, show_colnames = FALSE,
         filename = "results/ATAC/QC/ATAC_spearman_heatmap.pdf", width = 10, height = 8)

#Make a PCA plot of log2 TPM values
pca = performPCA(atac_tpm_log2, atac_data$design)
pca_plot = ggplot(pca$pca_matrix, aes(x = PC1, y = PC2, color = condition_name, label = sample_id)) + 
  geom_point()
ggsave("results/ATAC/QC/ATAC_tpm_pca.pdf", pca_plot, width = 5.5, height = 4)
