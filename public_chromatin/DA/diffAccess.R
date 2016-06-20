library("dplyr")
library("DESeq2")
load_all("macrophage-chromatin/housekeeping/")

#Construct design matrix
sample_names = read.table("macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
sample_names = sample_names[grepl("H3K27Ac",sample_names)]
design = constructDesignMatrix_Ivashkiv(sample_names)
design$condition = factor(design$condition, levels = c("naive","IFNg", "LPS", "IFNg_LPS"))
rownames(design) = design$sample_id

#Import peak coordinates
h3k27ac_peaks = import.gff3("annotations/H3K27Ac_joint_peaks.gff3")

#Import peak counts
counts = readRDS("results/Ivashkiv/H3K27Ac_combined_counts.rds")
length_df = dplyr::select(counts, gene_id, length)
rownames(counts) = counts$gene_id
diff_counts = counts[,rownames(design)]

dge = DGEList(counts = diff_counts)
dge = calcNormFactors(dge, method = "TMM")

#Apply voom transformation
design_matrix = model.matrix(~condition, design)
v <- voom(dge,design_matrix,plot=TRUE)
fit = lmFit(v, design_matrix)
fit <- eBayes(fit)

#Extract DE peaks
ifng_hits = topTable(fit,coef=c(2), lfc = 1, p.value = 0.1, number = 1e5) %>% tidyTopTable()
lps_hits = topTable(fit,coef=c(3), lfc = 1, p.value = 0.1, number = 1e5) %>% tidyTopTable()
ifng_lps_hits = topTable(fit,coef=c(4), lfc = 1, p.value = 0.1, number = 1e5) %>% tidyTopTable()
variable_peaks = unique(c(ifng_hits$gene_id, ifng_lps_hits$gene_id, lps_hits$gene_id))

#Calculate TPM of the expression level
length_df = dplyr::select(counts, gene_id, length)
tpm_matrix = log(calculateTPM(diff_counts, length_df) + 0.01, 2)

#Calculate mean peak height per condition
mean_tpm = calculateMean(tpm_matrix, as.data.frame(design), "condition")
mean_tpm = mean_tpm[variable_peaks,]

#Standardise effect sizes
cqn_set = ExpressionSet(as.matrix(mean_tpm))
cqn_set_std = standardise(cqn_set)

#Cluster the expression data
clusters = mfuzz(cqn_set_std, c = 5, m = 1.5, iter = 1000)
cluster_cores = acore(cqn_set_std, clusters, min.acore = 0.7)
names(cluster_cores) = c(1:length(cluster_cores))
pheatmap(clusters$centers)

#Make a df of clustered data
cluster_df = plyr::ldply(cluster_cores, .id = "cluster_id") %>%
  dplyr::rename(gene_id = NAME)
exp_df = exprs(cqn_set_std) %>% as.data.frame() %>% 
  dplyr::mutate(gene_id = rownames(exprs(cqn_set_std))) %>% tidyr::gather("condition_name", "expression", 1:4) %>%
  dplyr::mutate(condition_name = factor(as.character(condition_name), levels = c("naive","IFNg","LPS","IFNg_LPS")))
cluster_exp = dplyr::left_join(cluster_df, exp_df, by = "gene_id")

#Calculate means and order clusters
cluster_means = dplyr::group_by(cluster_exp, cluster_id, condition_name) %>% dplyr::summarise(expression = mean(expression))
cluster_order = tidyr::spread(cluster_means, condition_name, expression) %>% dplyr::arrange(naive) %>% 
  dplyr::mutate(naive_bin = ifelse(naive < 0, 0, 1)) %>% 
  dplyr::mutate(IFNg_bin = ifelse(IFNg < 0, 0, 1)) %>%
  dplyr::mutate(LPS_bin = ifelse(IFNg < LPS, 0, 1)) %>%
  dplyr::mutate(IFNg_LPS_bin = ifelse(IFNg_LPS < 0, 0, 1)) %>%
  arrange(naive_bin, IFNg_bin, LPS_bin) %>%
  dplyr::mutate(new_cluster_id = c(1:5)) %>%
  dplyr::select(cluster_id, new_cluster_id)

cluster_plot_data = dplyr::left_join(cluster_exp, cluster_order, by = "cluster_id") %>% 
  dplyr::group_by(new_cluster_id)

ggplot(cluster_plot_data, aes(x = condition_name, y = expression)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_boxplot()

#Make a heatmap of the H3K27Ac clusters
h3k27ac_heatmap = ggplot(cluster_plot_data %>% dplyr::sample_frac(0.5), aes(x = condition_name, y = gene_id, fill = expression)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Expression", midpoint = 0) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave("results/Ivashkiv/DA/H3K27Ac_heatmap.png", h3k27ac_heatmap, width  =5, height = 7)

#Give names to clusters
peak_coords = dplyr::select(as.data.frame(h3k27ac_peaks), seqnames, start, end, strand, type, gene_id)
cluster_names = data_frame(new_cluster_id = c(1:5), name = c("IFNg_LPS_up", "LPS_up","IFNg_up", "IFNg_down", "LPS_down"))
cluster_df = dplyr::select(cluster_plot_data, gene_id, new_cluster_id) %>% unique() %>% 
  dplyr::arrange(new_cluster_id) %>% 
  dplyr::left_join(cluster_names, by = "new_cluster_id") %>%
  dplyr::ungroup() %>%
  dplyr::select(-new_cluster_id) %>%
  dplyr::left_join(peak_coords, by = "gene_id") %>%
  dataFrameToGRanges()
export.bed(cluster_df, "results/Ivashkiv/DA/H3K27Ac_clustered_peaks.bed")

