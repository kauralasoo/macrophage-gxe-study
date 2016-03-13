library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")
library("readr")
library("limma")
library("edgeR")
library("Biobase")
library("Mfuzz")
library("ggplot2")
library("pheatmap")

#Import atac data list
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))

#Apply TMM normalisation
dge = DGEList(counts = atac_list$counts)
dge = calcNormFactors(dge, method = "TMM")

#Apply voom transformation
design_matrix = model.matrix(~condition_name, atac_list$sample_metadata)
v <- voom(dge,design_matrix,plot=TRUE)
fit = lmFit(v, design_matrix)
fit <- eBayes(fit)

#Extract signifcant peaks
ifng_hits = topTable(fit,coef=c(2), lfc = 2, p.value = 0.01, number = 1e5) %>% tidyTopTable()
sl1344_hits = topTable(fit,coef=c(3), lfc = 2, p.value = 0.01, number = 1e5) %>% tidyTopTable()
ifng_sll1344_hits = topTable(fit,coef=c(4), lfc = 2, p.value = 0.01, number = 1e5) %>% tidyTopTable()
variable_peaks = unique(c(ifng_hits$gene_id, ifng_sll1344_hits$gene_id, sl1344_hits$gene_id))

#Calculate mean peak height per condition
mean_tpm = calculateMean(atac_list$tpm, as.data.frame(atac_list$sample_metadata), "condition_name")
mean_tpm = log(mean_tpm + 0.01,2)
mean_tpm = mean_tpm[variable_peaks,]

#Standardise effect sizes
cqn_set = ExpressionSet(as.matrix(mean_tpm))
cqn_set_std = standardise(cqn_set)

#Cluster the expression data
clusters = mfuzz(cqn_set_std, c = 6, m = 1.5, iter = 1000)
cluster_cores = acore(cqn_set_std, clusters, min.acore = 0.7)
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
  dplyr::mutate(new_cluster_id = c(1:6)) %>%
  dplyr::select(cluster_id, new_cluster_id)

cluster_plot_data = dplyr::left_join(cluster_exp, cluster_order, by = "cluster_id") %>% 
  dplyr::group_by(new_cluster_id)

#Make a heatmap
chromatin_clusters_heatmap = ggplot(cluster_plot_data %>% dplyr::sample_frac(0.5), aes(x = condition_name, y = gene_id, fill = expression)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Expression", midpoint = 0) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave("results/ATAC/DA/ATAC_DA_clusters.png",chromatin_clusters_heatmap, width = 5, height = 7)

ggplot(cluster_plot_data, aes(x = condition_name, y = expression)) + 
  facet_grid(new_cluster_id ~ .,  scales = "free_y", space = "free_y") + stat_summary(fun.y = mean, geom="point")

#Save clusters to disk
final_clusters = dplyr::select(cluster_plot_data, gene_id, MEM.SHIP, new_cluster_id) %>% unique()
saveRDS(final_clusters, "results/ATAC/DA/peak_clusters.rds")
final_clusters = readRDS("results/ATAC/DA/peak_clusters.rds")

#Add coordinates to each cluster
cluster_names = data_frame(new_cluster_id = c(1:6), name = c("IFNg_SL1344_up", "SL1344_up", "inflammatory_up", "IFNg_up", "IFNg_down", "SL1344_down"))

#Export all clusters as a single bed file
cluster_ranges = dplyr::left_join(final_clusters, cluster_names, by = "new_cluster_id") %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(atac_list$gene_metadata) %>% 
  dplyr::transmute(gene_id, name, seqnames = chr, start, end, strand)
export.bed(dataFrameToGRanges(cluster_ranges), "results/ATAC/DA/ATAC_clustered_peaks.bed")

#Export ChIP peaks
chip_peaks = readRDS("results/ATAC/DA/Chip_peak_lists.rds")
export.bed(chip_peaks$Ivashkiv$STAT1_rep1_B, "results/ATAC/DA/ChIP_peak_lists/STAT1_rep1_B.bed")
export.bed(chip_peaks$Ivashkiv$STAT1_rep1_D, "results/ATAC/DA/ChIP_peak_lists/STAT1_rep1_D.bed")
export.bed(chip_peaks$Ivashkiv$IRF1_B, "results/ATAC/DA/ChIP_peak_lists/IRF1_B.bed")
export.bed(chip_peaks$Rehli$MAC_PU1, "results/ATAC/DA/ChIP_peak_lists/MAC_PU1.bed")
export.bed(chip_peaks$Rehli$MAC_CEBPbeta, "results/ATAC/DA/ChIP_peak_lists/MAC_CEBPbeta.bed")


gat-run.py --ignore-segment-tracks --segments=ChIP_peak_lists/STAT1_rep1_D.bed --annotations=ATAC_peak_lists/synergisitic.bed --workspace=../../../../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=10 --log=gat.log > gat_salmonella.out
gat-compare


#Perform motif enrichment
db <- CisBP.extdata("Homo_sapiens")
tfs <- tfbs.createFromCisBP(db)

#Agnes clustering
tfs1 <- tfbs.clusterMotifs(tfs, method="agnes", group.k=100, pdf.heatmap="agnes.hm.pdf")
#AP clustering
tfs2 <- tfbs.clusterMotifs(tfs, method="apcluster", pdf.heatmap= "apcluster.hm.pdf")

# Draw motif logos with one group of TF per page
tfbs.drawLogosForClusters(tfs1, file.pdf="agnes.logos.pdf");
tfbs.drawLogosForClusters(tfs2, file.pdf="apcluster.logos.pdf")


