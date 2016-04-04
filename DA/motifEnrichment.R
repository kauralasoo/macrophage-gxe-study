library("plyr")
library("dplyr")
library("readr")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

#Import ATAC data and clusters
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))
final_clusters = readRDS("results/ATAC/DA/peak_clusters.rds")

#Import FIMO motif matches
fimo_hits = readr::read_delim("results/ATAC/FIMO_CISBP_results.long.txt", delim = "\t", col_types = c("cciicddcc"), 
                              col_names = c("motif_id","seq_name","start","end","strand","score","p_value","dummy","matched_seq"), skip = 1)
fimo_hits_clean = tidyr::separate(fimo_hits, seq_name, c("prefix","gene_id"), sep = "=") 

#Filter motifs by TF expression
combined_expression_data = readRDS("../macrophage-gxe-study/results/SL1344/combined_expression_data.rds")
mean_tpm = calculateMean(combined_expression_data$tpm, as.data.frame(combined_expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_tpm, 1, max) > 1))
expressed_motifs = dplyr::filter(unique_motifs, gene_id %in% expressed_genes)

#Make mean_tpm df
mean_df = dplyr::mutate(mean_tpm, gene_id = rownames(mean_tpm)) %>%
  dplyr::left_join(dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name), by = "gene_id")

#Import motif metadata
TF_information = readr::read_tsv("~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/TF_Information.txt")
colnames(TF_information)[6] = "gene_id"
unique_motifs = dplyr::select(TF_information, Motif_ID, gene_id, TF_Name) %>% dplyr::filter(Motif_ID != ".") %>%
  dplyr::rename(motif_id = Motif_ID, tf_name = TF_Name)


#Calculate motif enrichments in each cluster
cluster_list = plyr::dlply(final_clusters,.(name))
cluster_enrichments = purrr::map(cluster_list, ~fimoRelativeEnrichment(.,NULL, fimo_hits_clean, atac_list$gene_metadata))
filtered_enrichments = purrr::map(cluster_enrichments, ~dplyr::left_join(.,unique_motifs, by = "motif_id") %>%
                                   dplyr::semi_join(expressed_motifs, by = "motif_id") %>%
                                   dplyr::arrange(-enrichment) )
motif_enrichment_df = ldply(filtered_enrichments, .id = "cluster_name") %>% tbl_df()

#Make a heatmap of motif enrichments
interesting_tfs = c("FOSB","BATF3","POU2F1","NFKB1","IRF1","STAT1", "MAFB", "MEF2A")
selected_enrichments = dplyr::filter(motif_enrichment_df, tf_name %in% interesting_tfs) %>%
  dplyr::mutate(l2fold = log(enrichment, 2)) %>%
  dplyr::mutate(cluster_name = factor(as.character(cluster_name), 
          levels =rev(c("SL1344_up_1","SL1344_up_2" ,"IFNg_SL1344_up", "IFNg_up_1", "IFNg_up_2", "IFNg_down", "SL1344_down")))) %>%
  dplyr::mutate(tf_name = factor(tf_name, levels = interesting_tfs))

#Make a plot of motif enrichments
motif_enrichment_plot = ggplot2::ggplot(selected_enrichments, aes(x = tf_name, y = cluster_name, fill = l2fold, label = round(l2fold,1))) + 
    geom_tile() +
    scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", name = "Log2 enrichment", midpoint = 0) +
    geom_text() +
    xlab("TF motif") + 
    ylab("ATAC-seq peak cluster") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
ggsave("results/ATAC/DA/Motif_enrichment_in_clusters.pdf", plot = motif_enrichment_plot, width = 8, height = 8)






