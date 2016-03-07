library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
library("DESeq2")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data.rds")
combined_expression_data_filtered = readRDS("results/SL1344/combined_expression_data_covariates.rds")
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)


design = dplyr::filter(combined_expression_data_filtered$sample_metadata, donor %in% c("yemz","ieki","bima","aipt","zuta")) %>%
  dplyr::arrange(condition) %>% as.data.frame()
rownames(design) = design$sample_id

count_matrix = combined_expression_data_filtered$counts[,design$sample_id]

cqn_set = ExpressionSet(as.matrix(mean_cqn))
cqn_set_std = standardise(cqn_set)

c1 = mfuzz(cqn_set_std, c = 6, m = 2)
mfuzz.plot(cqn_set_std, cl = c1, mfrow=c(4,4))
b = acore(cqn_set_std, c1, min.acore = 0.7)

#Use DESeq to identify genes that vary between conditions
dds = DESeq2::DESeqDataSetFromMatrix(count_matrix, design, ~condition_name) 
dds = DESeq(dds, test = "LRT", reduced = ~ 1)

ifng_genes = results(dds, contrast=c("condition_name","naive","IFNg")) %>% 
  tidyDESeq(gene_name_map) %>% dplyr::filter(padj < 0.01, abs(log2FoldChange) > 1)
sl1344_genes = results(dds, contrast=c("condition_name","naive","SL1344")) %>% 
  tidyDESeq(gene_name_map) %>% dplyr::filter(padj < 0.01, abs(log2FoldChange) > 1)
ifng_sl1344_genes = results(dds, contrast=c("condition_name","naive","IFNg_SL1344")) %>% 
  tidyDESeq(gene_name_map) %>% dplyr::filter(padj < 0.01, abs(log2FoldChange) > 1)
variable_genes = c(ifng_genes$gene_id, sl1344_genes$gene_id, ifng_sl1344_genes$gene_id) %>% unique()

#Calculate mean expression
cqn_matrix = combined_expression_data_filtered$cqn[,design$sample_id]
mean_cqn = calculateMean(cqn_matrix, design, "condition_name")
mean_cqn = mean_cqn[variable_genes,]

#Convert to Expression Set and standardise
cqn_set = ExpressionSet(as.matrix(mean_cqn))
cqn_set_std = standardise(cqn_set)

clusters = mfuzz(cqn_set_std, c = 9, m = 1.5, iter = 1000)
cluster_cores = acore(cqn_set_std, clusters, min.acore = 0.7)
names(cluster_cores) = c(1:length(cluster_cores))
pheatmap(clusters$centers)

cluster_df = ldply(cluster_cores, .id = "cluster_id") %>%
  dplyr::rename(gene_id = NAME)
exp_df = exprs(cqn_set_std) %>% as.data.frame() %>% 
  dplyr::mutate(gene_id = rownames(mean_cqn)) %>% tidyr::gather("condition_name", "expression", 1:4) %>%
  dplyr::mutate(condition_name = factor(as.character(condition_name), levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
cluster_exp = dplyr::left_join(cluster_df, exp_df, by = "gene_id")

#Calculate means for eacg cluster
cluster_means = dplyr::group_by(cluster_exp, cluster_id, condition_name) %>% summarise(expression = mean(expression))
cluster_order = tidyr::spread(cluster_means, condition_name, expression) %>% arrange(naive) %>% 
  dplyr::mutate(naive_bin = ifelse(naive < 0, 0, 1)) %>% 
  dplyr::mutate(IFNg_bin = ifelse(IFNg < 0, 0, 1)) %>%
  dplyr::mutate(SL1344_bin = ifelse(IFNg < SL1344, 0, 1)) %>%
  dplyr::mutate(IFNg_SL1344_bin = ifelse(IFNg_SL1344 < 0, 0, 1)) %>%
  arrange(naive_bin, IFNg_bin, SL1344_bin) %>%
  dplyr::mutate(new_cluster_id = c(1:9)) %>%
  dplyr::select(cluster_id, new_cluster_id)

cluster_exp = dplyr::left_join(cluster_exp, cluster_order)

diff_exp_heatmap = ggplot(cluster_exp, aes(x = condition_name, y = gene_id, fill = expression)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Expression", midpoint = 0) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave("results/SL1344/SL1344_DE_clusters.pdf",diff_exp_heatmap, width = 5, height = 7)


