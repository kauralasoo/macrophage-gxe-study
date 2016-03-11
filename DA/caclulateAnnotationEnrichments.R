#Run GAT
#gat-run.py --segments=ATAC/DA/ATAC_clustered_peaks.bed --annotations=Ivashkiv/DA/H3K27Ac_clustered_peaks.bed --workspace=../../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=gat.log --with-segment-tracks > gat.out

#Import GAT results
gat_results = read.table("results/gat.out", header = TRUE, stringsAsFactors = FALSE)
gat_results$track = factor(gat_results$track, 
        levels = c("IFNg_up", "SL1344_up", "inflammatory_up", "IFNg_SL1344_up", "IFNg_down", "SL1344_down"))
gat_results$annotation = factor(gat_results$annotation, 
        levels = c("IFNg_up", "LPS_up", "IFNg_LPS_up", "IFNg_down", "LPS_down"))

peak_enrichments = ggplot(gat_results, aes(x = track, y = annotation, fill = l2fold, label = round(l2fold,1))) + 
  geom_tile() +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Log2 enrichment", midpoint = 0) +
  geom_text() +
  xlab("ATAC-Seq peaks (this study)") + 
  ylab("H3K27Ac peaks (Qiao et al)") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

ggsave("results/ATAC/DA/MDM_vs_IPSDM_concordance.pdf", peak_enrichments, width = 10, height = 7)
