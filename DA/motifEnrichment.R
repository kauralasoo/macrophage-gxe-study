library("dplyr")
library("readr")

#Import ATAC data and clusters
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))
final_clusters = readRDS("results/ATAC/DA/peak_clusters.rds")

#Import FIMO motif matches
fimo_hits = readr::read_delim("results/ATAC/FIMO_CISBP_results.long.txt", delim = "\t", col_types = c("cciicddcc"), 
                              col_names = c("motif_id","seq_name","start","end","strand","score","p_value","dummy","matched_seq"), skip = 1)

#Count the number of matches per peak
fimo_hits_clean = tidyr::separate(fimo_hits,seq_name, c("prefix","gene_id"), sep = "=") 
width_matrix = dplyr::mutate(atac_list$gene_metadata, width = end - start)

#Import motif metadata
TF_information = readr::read_tsv("~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/TF_Information.txt")
colnames(TF_information)[6] = "gene_id"
unique_motifs = dplyr::select(TF_information, Motif_ID, gene_id, TF_Name) %>% dplyr::filter(Motif_ID != ".") %>%
  dplyr::rename(motif_id = Motif_ID, tf_name = TF_Name)


#Calculate total number of matches for each motif
total_matches = dplyr::group_by(fimo_hits_clean,motif_id) %>%
  dplyr::summarise(n_matches = length(motif_id)) %>%
  dplyr::mutate(total_length = sum(width_matrix$width))

sl1344_up_peaks = dplyr::filter(final_clusters, name == "IFNg_up")
sl1344_lengths = sum(dplyr::filter(width_matrix, gene_id %in% sl1344_up_peaks$gene_id)$width)
sl1344_matches = dplyr::filter(fimo_hits_clean, gene_id %in% sl1344_up_peaks$gene_id) %>% 
  dplyr::group_by(motif_id) %>% 
  summarise(n_cluster_matches = length(motif_id)) %>%
  dplyr::mutate(cluster_length = sl1344_lengths)

#Motif_enrichment
a = dplyr::left_join(total_matches, sl1344_matches, by = "motif_id") %>% 
  dplyr::mutate(enrichment = (n_cluster_matches/cluster_length)/(n_matches/total_length)) %>% 
  arrange(-enrichment) %>% 
  dplyr::group_by(motif_id) %>%
  dplyr::mutate(p_hyper = phyper(n_cluster_matches - 1, n_matches, total_length - n_matches, cluster_length, lower.tail = FALSE)) %>%
  dplyr::mutate(p_fdr = p.adjust(p_hyper, "fdr")) %>%
  dplyr::left_join(unique_motifs, by = "motif_id")
