library("dplyr")
library("pheatmap")
library("devtools")
library("ggplot2")
load_all("macrophage-chromatin/housekeeping/")
load_all("../seqUtils/")

#Import peak counts
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
atac_data$sample_metadata$condition_name = factor(atac_data$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))

#Construct design matrix
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
design = constructDesignMatrix_ATAC(sample_names)

#Calculate TPM values
atac_tpm_log2 = log(atac_data$tpm + 0.1, 2)

#Make heatmap of TPM counts
pheatmap(spearman_corr[filtered_metadata$sample_id,filtered_metadata$sample_id], width = 10, height = 10)

#Make a heatmap of CQN normalized accesibility data
spearman_corr = cor(atac_tpm_log2,method = "spearman")
meta = as.data.frame(dplyr::select(atac_data$sample_metadata, condition_name))
rownames(meta) = atac_data$sample_metadata$sample_id
pheatmap(spearman_corr, annotation_row = meta, 
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         filename = "results/ATAC/QC/ATAC_spearman_heatmap.pdf", width = 10, height = 8)

#Make a PCA plot of log2 TPM values
pca = performPCA(atac_tpm_log2, atac_data$sample_metadata)
var_exp = pca$var_exp * 100 
xlabel = paste0("PC1 (", round(var_exp[1]), "%)")
ylabel = paste0("PC2 (", round(var_exp[2]), "%)")
pca_plot = ggplot(pca$pca_matrix, aes(x = PC1, y = -PC2, color = condition_name, label = sample_id)) + 
  geom_point() + 
  xlab(xlabel) +
  ylab(ylabel) + 
  scale_color_manual(values = conditionPalette(), name = "condition") + 
  theme_light() + 
  theme(legend.key = element_blank())
ggsave("results/ATAC/DA/PCA_of_chromatin_accessibility.pdf", pca_plot, width = 5.5, height = 4)

#Perform PCA within each condition seprately
IFNg_design = dplyr::filter(atac_data$sample_metadata, condition_name == "IFNg")
IFNg_tpm = atac_tpm_log2[,IFNg_design$sample_id]
IFNg_pca = performPCA(IFNg_tpm, IFNg_design)
ggplot(IFNg_pca$pca_matrix, aes(x = PC3, y = PC4, color = condition_name, label = sample_id)) + 
  geom_point() + geom_text()

IFNg_design = dplyr::filter(atac_data$sample_metadata, condition_name == "IFNg_SL1344")
IFNg_tpm = atac_tpm_log2[,IFNg_design$sample_id]
IFNg_pca = performPCA(IFNg_tpm, IFNg_design)
ggplot(IFNg_pca$pca_matrix, aes(x = PC1, y = PC2, color = condition_name, label = sample_id)) + 
  geom_point() + geom_text()
