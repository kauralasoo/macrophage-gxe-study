library("ggplot2")
library("devtools")
load_all("../seqUtils/")

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
gat_results = read.table("results/public_chromatin/annotation_overlaps/H3K27Ac_overlap.gat.txt", header = TRUE, stringsAsFactors = FALSE)
gat_results$track = factor(gat_results$track, 
        levels = rev(c("Salmonella","Both" ,"IFNg", "Decreased")))
gat_results$annotation = factor(gat_results$annotation, 
        levels = rev(c("LPS", "Both", "IFNg", "Decreased")))

peak_enrichments = plotEnrichment(gat_results, ylabel = "ATAC-Seq peaks (this study)", xlabel = "H3K27Ac peaks (Qiao et al)")
ggsave("figures/supplementary/concordance_MDM_vs_IPSDM.pdf", peak_enrichments, width = 4.5, height = 4.5)


#CIITA-RFX5
ciita_enrich = read.table("results/public_chromatin/annotation_overlaps/CIITA-RFX5_overlap.gat.txt", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(track, annotation, l2fold) %>%
  dplyr::filter(annotation %in% c("CIITA_IFNg", "RFX5_IFNg"))

#Make a plot
ciita_plot = ggplot(ciita_enrich, aes(x = l2fold, y = annotation)) + 
  geom_point() + 
  facet_grid(track~.) + 
  theme_light() + 
  geom_vline(xintercept = 0) + 
  xlab("Log2 fold-enrichment") +
  ylab("Monocyte ChIP-seq data")
ggsave("figures/supplementary/concordance_CIITA-RFX5.pdf", ciita_plot, width = 3.5, height = 4.5)

#PU.1 and CEBPb from two different papers
naive_enrich = read.table("results/ATAC/ChIP_enrichment/gat_output/naive_naive_overlap.gat.txt", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(track, annotation, l2fold)
plotEnrichment(naive_enrich, xlabel = "Monocyte ChIP-Seq", ylabel = "ATAC-Seq peaks")

