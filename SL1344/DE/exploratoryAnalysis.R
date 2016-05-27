library("DESeq2")
library("dplyr")
library("devtools")
library("tidyr")
library("ggplot2")
library("gplots")
library("cqn")
library("gProfileR")
load_all("../seqUtils/")
library("pheatmap")

#Load expression data from disk
eqtl_data_list = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Make heatmap
#Extract metadata for the heatmap command
meta = as.data.frame(dplyr::select(eqtl_data_list$sample_metadata, condition_name))
rownames(meta) = eqtl_data_list$sample_metadata$sample_id

pheatmap(cor(eqtl_data_list$cqn, method = "spearman"), annotation_row = meta, 
         show_rownames = TRUE, show_colnames = FALSE,
         filename = "results/SL1344/diffExp/sample_heatmap.pdf", width = 20, height = 20)
pheatmap(cor(eqtl_data_list$cqn, method = "spearman"), annotation_row = meta, 
         show_rownames = FALSE, show_colnames = FALSE)

#Perform PCA analysis
pca_list = performPCA(eqtl_data_list$cqn, eqtl_data_list$sample_metadata)
var_exp = pca_list$var_exp * 100 
xlabel = paste0("PC1 (", round(var_exp[1]), "%)")
ylabel = paste0("PC2 (", round(var_exp[2]), "%)")
pca_list$pca_matrix = dplyr::mutate(pca_list$pca_matrix, 
    condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344", "IFNg_SL1344")))
pca_plot = ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition_name)) + 
  geom_point() +
  xlab(xlabel) + 
  ylab(ylabel) +
  scale_color_manual(values = conditionPalette(), name = "condition") + 
  theme_light() + 
  theme(legend.key = element_blank())
ggsave("results/SL1344/DE/PCA_of_gene_expression.pdf", plot = pca_plot, width = 5.5, height = 4)

#Make plot of IFNb
ifnb_plot = plotGene("ENSG00000171855",expression_list$exprs_cqn, design_matrix, expression_list$gene_metadata)
ggsave("results/SL1344/diffExp/IFNB1.pdf", plot = ifnb_plot, width = 5, height = 5)

#Make some plots for subho
subho_genes = read.table("../LPS/subho_genes.txt", stringsAsFactors = FALSE)[,1]
plots = lapply(as.list(subho_genes), plotGene, expression_list$exprs_cqn, design_matrix, expression_list$gene_metadata)
savePlots(plots, "subho_plots/", 6, 6)
