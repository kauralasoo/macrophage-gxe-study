library("DESeq2")
library("dplyr")
library("devtools")
library("tidyr")
library("ggplot2")
library("gplots")
library("cqn")
library("pheatmap")
load_all("../seqUtils/")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")


#Load expression data from disk
eqtl_data_list = readRDS("results/SL1344/combined_expression_data_covariates.rds")
eqtl_data_list$sample_metadata$condition_name = factor(eqtl_data_list$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))

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
pca_data = dplyr::mutate(pca_list$pca_matrix, 
    condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344", "IFNg_SL1344"))) %>%
  dplyr::left_join(figureNames(), by = "condition_name")
pca_plot = ggplot(pca_data, aes(x = PC1, y = PC2, color = figure_name)) + 
  geom_point() +
  xlab(xlabel) + 
  ylab(ylabel) +
  scale_color_manual(values = conditionPalette(), name = "Condition") + 
  theme_light() + 
  theme(legend.key = element_blank())
ggsave("figures/supplementary/PCA_of_gene_expression.pdf", plot = pca_plot, width = 4.5, height = 3.5)

#Look at the effect of GC bias
pca_list = performPCA(log(eqtl_data_list$tpm + 0.01,2), eqtl_data_list$sample_metadata)
pca_list2 = performPCA(eqtl_data_list$cqn, eqtl_data_list$sample_metadata)

pca_plot = ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC5, color = rna_auto)) + 
  geom_point() +
  theme_light() + 
  theme(legend.key = element_blank())
ggplot(pca_list2$pca_matrix, aes(x = PC1, y = PC5, color = rna_auto)) + 
  geom_point() +
  theme_light() + 
  theme(legend.key = element_blank())

#Make plot of IFNb
ifnb_plot = plotGene("ENSG00000179583",eqtl_data_list$cqn, eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)
ggsave("results/SL1344/diffExp/IFNB1.pdf", plot = ifnb_plot, width = 5, height = 5)

#Make some plots for subho
subho_genes = read.table("../LPS/subho_genes.txt", stringsAsFactors = FALSE)[,1]
plots = lapply(as.list(subho_genes), plotGene, expression_list$exprs_cqn, design_matrix, expression_list$gene_metadata)
savePlots(plots, "subho_plots/", 6, 6)

#Make a plot of the CIITA gene
ciita_data = extractGeneData("ENSG00000179583",eqtl_data_list$cqn, eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata) %>%
  dplyr::left_join(friendlyNames())

ciita_plot = ggplot(ciita_data, aes(x = friendly_name, y = expression, fill = friendly_name)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = .2)) + 
  ylab(expression(paste(Log[2], " normalised expression"))) + 
  ggtitle(ciita_data$gene_name[1]) +
  scale_fill_manual(values = conditionPalette()) + 
  theme_light() + 
  theme(legend.position="none", axis.title.x = element_blank())
ggsave("figures/supplementary/CIITA_expression.pdf", ciita_plot, width = 4, height = 4)

#Explore total read counts
eqtl_data_list_full = readRDS("results/SL1344/combined_expression_data.rds")
fc_total = colSums(eqtl_data_list_full$counts) %>% tidyVector(value_id = "fragment_count")
fc_df = dplyr::left_join(eqtl_data_list_full$sample_metadata, fc_total, by = "sample_id") %>%
  dplyr::mutate(library_type = ifelse(rna_auto, "automatic", "manual")) %>%
  dplyr::mutate(fragment_count = fragment_count/1e6)

#Make a histogram of library size
fragment_count_plot = ggplot(fc_df, aes(x = fragment_count, fill = library_type)) + 
  geom_histogram(position = "identity", alpha = .6) + 
  theme_light() + 
  theme(legend.key = element_blank()) + 
  xlab("Library size (millions of fragments)") +
  ylab("Sample count") +
  theme(legend.position = "top")
ggsave("figures/supplementary/rna_library_size_by_protocol.pdf", 
       plot = fragment_count_plot, width = 4, height = 4)

#Explore GC bias
#Explore the GC bias across samples
expression_cqn = cqn(counts = eqtl_data_list_full$counts[eqtl_data_list_full$gene_metadata$gene_id,], 
                     x = eqtl_data_list_full$gene_metadata$percentage_gc_content, 
                     lengths = eqtl_data_list_full$gene_metadata$length, verbose = TRUE)

is_rna_auto = dplyr::select(eqtl_data_list_full$sample_metadata, sample_id, rna_auto)

gc_grid = dplyr::mutate(as.data.frame(expression_cqn$func1), 
                        grid_point = expression_cqn$grid1) %>%
  tidyr::gather(sample_id, gc_effect, aipt_A:zuta_D) %>%
  dplyr::left_join(is_rna_auto, by = "sample_id") %>%
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(gc_effect = gc_effect - mean(gc_effect)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(library_type = ifelse(rna_auto, "automatic", "manual")) %>%
  tbl_df()

gc_bias = ggplot(gc_grid, aes(x = grid_point, y = gc_effect, group = library_type, color = library_type)) +
  stat_smooth() +
  theme_light() + 
  theme(legend.key = element_blank()) +
  xlab("GC content (%)") +
  ylab("Relative bias") +
  theme(legend.position = "top")
ggsave("figures/supplementary/rna_gc_bias_by_protocol.pdf", 
       plot = gc_bias, width = 4, height = 4)

