library("ggplot2")
#Run GAT
#gat-run.py --segments=ATAC/DA/ATAC_clustered_peaks.bed --annotations=Ivashkiv/DA/H3K27Ac_clustered_peaks.bed --workspace=../../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=gat.log --with-segment-tracks > gat.out

plotEnrichment <- function(df, xlabel = "", ylabel = ""){
  peak_enrichments = ggplot2::ggplot(df, aes(x = annotation, y = track, fill = l2fold, label = round(l2fold,1))) + 
    geom_tile() +
    scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Log2 enrichment", midpoint = 0) +
    geom_text() +
    xlab(xlabel) + 
    ylab(ylabel) + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
  return(peak_enrichments)
}

#Import GAT results
gat_results = read.table("results/annotation_overlaps/H3K27Ac_overlap.gat.txt", header = TRUE, stringsAsFactors = FALSE)
gat_results$track = factor(gat_results$track, 
        levels = rev(c("SL1344_up_1","SL1344_up_2" ,"IFNg_SL1344_up", "IFNg_up_1", "IFNg_up_2", "IFNg_down", "SL1344_down")))
gat_results$annotation = factor(gat_results$annotation, 
        levels = rev(c("LPS_up", "IFNg_LPS_up", "IFNg_up", "IFNg_down", "LPS_down")))

peak_enrichments = plotEnrichment(gat_results, ylabel = "ATAC-Seq peaks (this study)", xlabel = "H3K27Ac peaks (Qiao et al)")
ggsave("results/ATAC/DA/MDM_vs_IPSDM_concordance.pdf", peak_enrichments, width = 10, height = 7)


#Overlap enrichment with STAT1 peaks
gat-run.py --segments=results/ATAC/DA/ATAC_clustered_peaks.bed --annotations=results/Ivashkiv/DA/H3K27Ac_clustered_peaks.bed --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=results/annotation_overlaps/H3K27Ac_overlap.gat_log.txt --with-segment-tracks > results/annotation_overlaps/H3K27Ac_overlap.gat.txt
gat-run.py --segments=results/ATAC/DA/ATAC_clustered_peaks.bed --annotations=results/Ivashkiv/DA/STAT1_grouped_peaks.bed  --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=results/annotation_overlaps/STAT1_overlap.gat_log.txt --with-segment-tracks > results/annotation_overlaps/STAT1_overlap.gat.txt
gat-run.py --segments=results/ATAC/DA/ATAC_clustered_peaks.bed --annotations=results/Ivashkiv/peak_calls/IRF1_joint_peaks.bed --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=results/annotation_overlaps/IRF1_overlap.gat_log.txt --with-segment-tracks > results/annotation_overlaps/IRF1_overlap.gat.txt
gat-run.py --segments=results/ATAC/DA/ATAC_clustered_peaks.bed --annotations=results/Knight/CIITA-RFX5_joint_peaks.bed --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=results/annotation_overlaps/CIITA-RFX5_overlap.gat_log.txt --with-segment-tracks > results/annotation_overlaps/CIITA-RFX5_overlap.gat.txt

#Enrichments for PU.1, CEBPb and CTCF
gat-run.py --segments=annotations/ATAC_consensus_peaks.bed  --annotations=results/ATAC/ChIP_enrichment/naive_combined_peaks.bed --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=100 --log=results/ATAC/ChIP_enrichment/gat_output/ATAC_naive_overlap.gat_log.txt --with-segment-tracks > results/ATAC/ChIP_enrichment/gat_output/ATAC_naive_overlap.gat.txt
gat-run.py --segments=results/ATAC/ChIP_enrichment/naive_combined_peaks.bed  --annotations=results/ATAC/ChIP_enrichment/naive_combined_peaks.bed --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=100 --log=results/ATAC/ChIP_enrichment/gat_output/naive_naive_overlap.gat_log.txt --with-segment-tracks > results/ATAC/ChIP_enrichment/gat_output/naive_naive_overlap.gat.txt


#CIITA-RFX5
ciita_enrich = read.table("results/annotation_overlaps/CIITA-RFX5_overlap.gat.txt", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(track, annotation, l2fold)
ciita_enrich$track = factor(ciita_enrich$track, 
      levels = rev(c("SL1344_up_1","SL1344_up_2" ,"IFNg_SL1344_up", "IFNg_up_1", "IFNg_up_2", "IFNg_down", "SL1344_down")))
ciita_enrich_plot = plotEnrichment(ciita_enrich, xlabel = "Monocyte ChIP-Seq", ylabel = "ATAC-Seq peaks")
ggsave("results/ATAC/DA/CIITA-RFX5_concordance.pdf", ciita_enrich_plot, width = 7, height = 7)

#PU.1 and CEBPb from two different papers
naive_enrich = read.table("results/ATAC/ChIP_enrichment/gat_output/naive_naive_overlap.gat.txt", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(track, annotation, l2fold)
plotEnrichment(naive_enrich, xlabel = "Monocyte ChIP-Seq", ylabel = "ATAC-Seq peaks")

