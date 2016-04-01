library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
library("DESeq2")
library("Mfuzz")
library("pheatmap")
library("gProfileR")


#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data.rds")
combined_expression_data_filtered = readRDS("results/SL1344/combined_expression_data_covariates.rds")
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#### DESeq2 ####
#Use DESeq to identify genes that vary between conditions
design = combined_expression_data$sample_metadata %>% as.data.frame()
rownames(design) = design$sample_id
dds = DESeq2::DESeqDataSetFromMatrix(combined_expression_data$counts, design, ~condition_name) 
dds = DESeq(dds, test = "LRT", reduced = ~ 1)
saveRDS(dds, "results/SL1344/DE/DESeq2_condition_name_LRT_results.rds")

#Extract differentially expressed genes in each condition
ifng_genes = results(dds, contrast=c("condition_name","IFNg","naive")) %>% 
  tidyDESeq(gene_name_map) 
sl1344_genes = results(dds, contrast=c("condition_name","SL1344","naive")) %>% 
  tidyDESeq(gene_name_map) 
ifng_sl1344_genes = results(dds, contrast=c("condition_name","IFNg_SL1344","naive")) %>% 
  tidyDESeq(gene_name_map) 

#Construct a single matrix of fold-changed
log2FC_table = dplyr::transmute(ifng_genes, gene_id, IFNg_log2FC = log2FoldChange) %>%
  dplyr::left_join(dplyr::transmute(sl1344_genes, gene_id, SL1344_log2FC = log2FoldChange), by = "gene_id") %>%
  dplyr::left_join(dplyr::transmute(ifng_sl1344_genes, gene_id, IFNg_SL1344_log2FC = log2FoldChange), by = "gene_id") 

#Identify all differentially expressed genes
fc_genes = log2FC_table[which(apply(abs(log2FC_table[2:4]),1, max) > 1),]
variable_genes = intersect(dplyr::filter(ifng_genes, padj < 0.01)$gene_id, fc_genes$gene_id)


##### Clustering ####
#Claculate mean expression using TPM values
mean_tpm = calculateMean(combined_expression_data$tpm, design, "condition_name")
expressed_genes = names(which(apply(mean_tpm, 1, max) > 0.5))

#Calculate mean CQN and standardise
mean_cqn = calculateMean(combined_expression_data$cqn, design, "condition_name")
variable_cqn = mean_cqn[intersect(rownames(combined_expression_data_filtered$counts), variable_genes),]
cqn_set = ExpressionSet(as.matrix(variable_cqn))
cqn_set_std = standardise(cqn_set)

#Cluster the expression data
clusters = mfuzz(cqn_set_std, c = 9, m = 1.5, iter = 1000)
cluster_cores = acore(cqn_set_std, clusters, min.acore = 0.7)
names(cluster_cores) = c(1:length(cluster_cores))
pheatmap(clusters$centers)

cluster_df = ldply(cluster_cores, .id = "cluster_id") %>%
  dplyr::rename(gene_id = NAME)
exp_df = exprs(cqn_set_std) %>% as.data.frame() %>% 
  dplyr::mutate(gene_id = rownames(exprs(cqn_set_std))) %>% tidyr::gather("condition_name", "expression", 1:4) %>%
  dplyr::mutate(condition_name = factor(as.character(condition_name), levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
cluster_exp = dplyr::left_join(cluster_df, exp_df, by = "gene_id")

#Calculate means for each cluster
cluster_means = dplyr::group_by(cluster_exp, cluster_id, condition_name) %>% dplyr::summarise(expression = mean(expression))
cluster_order = tidyr::spread(cluster_means, condition_name, expression) %>% dplyr::arrange(naive) %>% 
  dplyr::mutate(naive_bin = ifelse(naive < 0, 0, 1)) %>% 
  dplyr::mutate(IFNg_bin = ifelse(IFNg < 0, 0, 1)) %>%
  dplyr::mutate(SL1344_bin = ifelse(IFNg < SL1344, 0, 1)) %>%
  dplyr::mutate(IFNg_SL1344_bin = ifelse(IFNg_SL1344 < 0, 0, 1)) %>%
  arrange(naive_bin, IFNg_bin, SL1344_bin) %>%
  dplyr::mutate(new_cluster_id = c(1:9)) %>%
  dplyr::select(cluster_id, new_cluster_id)

cluster_plot_data = dplyr::left_join(cluster_exp, cluster_order, by = "cluster_id") %>% 
  dplyr::group_by(new_cluster_id)

diff_exp_heatmap = ggplot(cluster_plot_data %>% dplyr::sample_frac(0.4), aes(x = condition_name, y = gene_id, fill = expression)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Expression", midpoint = 0) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave("results/SL1344/DE/DE_clusters.pdf",diff_exp_heatmap, width = 5, height = 7)
ggsave("results/SL1344/DE/DE_clusters.png",diff_exp_heatmap, width = 5, height = 7)

#Perform GO analysis
clusters = dplyr::select(cluster_plot_data, new_cluster_id, gene_id) %>% unique() %>%
  dplyr::left_join(log2FC_table, by = "gene_id") %>% dplyr::ungroup()

### Perform enrichment analysis with g:Profiler
#Salmonella both
cluster1 = (dplyr::filter(clusters, new_cluster_id == 1) %>% dplyr::arrange(-IFNg_SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster2 = (dplyr::filter(clusters, new_cluster_id == 2) %>% dplyr::arrange(-SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster3 = (dplyr::filter(clusters, new_cluster_id == 3) %>% dplyr::arrange(-IFNg_SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster4 = (dplyr::filter(clusters, new_cluster_id == 4) %>% dplyr::arrange(-IFNg_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster5 = (dplyr::filter(clusters, new_cluster_id == 5) %>% dplyr::arrange(-IFNg_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster6 = (dplyr::filter(clusters, new_cluster_id == 6) %>% dplyr::arrange(IFNg_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster7 = (dplyr::filter(clusters, new_cluster_id == 7) %>% dplyr::arrange(SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster8 = (dplyr::filter(clusters, new_cluster_id == 8) %>% dplyr::arrange(IFNg_SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster9 = (dplyr::filter(clusters, new_cluster_id == 9) %>% dplyr::arrange(IFNg_SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

gprofiler_results = list(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, cluster7, cluster8, cluster9)
saveRDS(gprofiler_results, "results/SL1344/DE/Mfuzz_clustering_go_enrichments.rds")