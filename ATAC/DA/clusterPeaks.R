library("devtools")
library("dplyr")
library("readr")
library("limma")
library("edgeR")
library("Biobase")
library("Mfuzz")
library("ggplot2")
library("pheatmap")
load_all("../seqUtils/")

#Import atac data list
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))

#Only keep samples with high quality data in all conditions
hq_donors = c("qolg","vass","kuxp","cicb", "febc","eiwy","oapg","nukw","hayt","bima","pamv","guss","eipl","iill","podx","pelm")
filtered_metadata = dplyr::filter(atac_list$sample_metadata, donor %in% hq_donors)

#Apply TMM normalisation
dge = DGEList(counts = atac_list$counts[,filtered_metadata$sample_id])
dge = calcNormFactors(dge, method = "TMM")

#Apply voom transformation
design_matrix = model.matrix(~condition_name, filtered_metadata)
v <- voom(dge,design_matrix,plot=TRUE)
fit = lmFit(v, design_matrix)
fit <- eBayes(fit)

#Extract signifcant peaks
ifng_hits = topTable(fit,coef=c(2), lfc = 2, p.value = 0.01, number = 1e5) %>% tidyTopTable()
sl1344_hits = topTable(fit,coef=c(3), lfc = 2, p.value = 0.01, number = 1e5) %>% tidyTopTable()
ifng_sll1344_hits = topTable(fit,coef=c(4), lfc = 2, p.value = 0.01, number = 1e5) %>% tidyTopTable()
variable_peaks = unique(c(ifng_hits$gene_id, ifng_sll1344_hits$gene_id, sl1344_hits$gene_id))

#Calculate mean peak height per condition
mean_cqn = calculateMean(atac_list$cqn[,filtered_metadata$sample_id], as.data.frame(filtered_metadata), "condition_name")
mean_cqn = mean_cqn[variable_peaks,]

#Standardise effect sizes
cqn_set = ExpressionSet(as.matrix(mean_cqn))
cqn_set_std = standardise(cqn_set)

#Cluster the expression data
set.seed(1)
clusters = mfuzz(cqn_set_std, c = 7, m = 1.5, iter = 1000)
cluster_cores = acore(cqn_set_std, clusters, min.acore = 0)
names(cluster_cores) = c(1:length(cluster_cores))
pheatmap(clusters$centers)

#Make a df of clustered data
cluster_df = plyr::ldply(cluster_cores, .id = "cluster_id") %>%
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
  dplyr::ungroup() %>%
  dplyr::mutate(new_cluster_id = c(1:7)) %>%
  dplyr::select(cluster_id, new_cluster_id)

cluster_plot_data = dplyr::left_join(cluster_exp, cluster_order, by = "cluster_id") %>% 
  dplyr::left_join(figureNames(), by = "condition_name") %>%
  dplyr::group_by(new_cluster_id)

#Make a heatmap
ylabel = paste(nrow(cluster_df), "accessible regions")
chromatin_clusters_heatmap = ggplot(cluster_plot_data %>% dplyr::sample_frac(0.4), 
                                    aes(x = figure_name, y = gene_id, fill = expression)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Accessibility", midpoint = 0) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x = element_blank(),legend.title = element_text(angle = 90)) + 
  theme(legend.position = "left") +
  ylab(ylabel) +
  theme(panel.spacing = unit(0.2, "lines"))
ggsave("figures/main_figures/ATAC_DA_clusters.png",chromatin_clusters_heatmap, width = 2.7, height = 4)

#Add names to each cluster
cluster_names = data_frame(new_cluster_id = c(1:7), name = c("Salmonella","Salmonella","Both", "IFNg", "IFNg", "Decreased", "Decreased"))

#Save clusters to disk
final_clusters = dplyr::select(cluster_plot_data, gene_id, MEM.SHIP, new_cluster_id) %>% unique() %>%
  dplyr::left_join(cluster_names, by = "new_cluster_id") %>% dplyr::ungroup()
saveRDS(final_clusters, "results/ATAC/DA/peak_clusters.rds")
final_clusters = readRDS("results/ATAC/DA/peak_clusters.rds") %>%
  dplyr::mutate(name = new_cluster_id)

#Export all clusters as a single bed file
cluster_ranges = final_clusters %>%
  dplyr::left_join(atac_list$gene_metadata) %>% 
  dplyr::transmute(gene_id, name, seqnames = chr, start, end, strand)

#Save clusters into a BED file for downstream analyis with GAT
rtracklayer::export.bed(dataFrameToGRanges(cluster_ranges), "results/ATAC/DA/ATAC_clustered_peaks.bed")
