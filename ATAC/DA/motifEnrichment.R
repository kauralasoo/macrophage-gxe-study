library("plyr")
library("dplyr")
library("readr")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")
library("seqLogo")

#Import ATAC data and clusters
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))
final_clusters = readRDS("results/ATAC/DA/peak_clusters.rds")

#Import FIMO motif matches
fimo_hits = readr::read_delim("results/ATAC/cisBP/FIMO_CISBP_results.long.txt.gz", delim = "\t", col_types = c("cciicddcc"), 
                              col_names = c("motif_id","seq_name","start","end","strand","score","p_value","dummy","matched_seq"), skip = 1)
fimo_hits_clean = tidyr::separate(fimo_hits, seq_name, c("prefix","gene_id"), sep = "=") 

#Import motif metadata
TF_information = readr::read_tsv("~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/TF_Information.txt")
colnames(TF_information)[6] = "gene_id"
unique_motifs = dplyr::select(TF_information, Motif_ID, gene_id, TF_Name) %>% dplyr::filter(Motif_ID != ".") %>%
  dplyr::rename(motif_id = Motif_ID, tf_name = TF_Name)

#Filter motifs by TF expression
combined_expression_data = readRDS("../macrophage-gxe-study/results/SL1344/combined_expression_data.rds")
mean_tpm = calculateMean(combined_expression_data$tpm, as.data.frame(combined_expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_tpm, 1, max) > 1))
expressed_motifs = dplyr::filter(unique_motifs, gene_id %in% expressed_genes)

#Make mean_tpm df
mean_df = dplyr::mutate(mean_tpm, gene_id = rownames(mean_tpm)) %>%
  dplyr::left_join(dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name), by = "gene_id") %>% tbl_df()

#Combine QTL clusters into larger groups to increase power
cluster_groups = data_frame(new_cluster_id = c(1:7), 
                            cluster_group_id = c("1-2","1-2","3","4-5","4-5","6-7","6-7"))
grouped_clusters = dplyr::left_join(final_clusters, cluster_groups, by = "new_cluster_id") %>%
  dplyr::transmute(gene_id, MEM.SHIP, new_cluster_id = cluster_group_id)

#Calculate motif enrichments in each cluster
cluster_list = plyr::dlply(grouped_clusters,.(new_cluster_id))
cluster_enrichments = purrr::map(cluster_list, ~fimoRelativeEnrichment(.,NULL, fimo_hits_clean, atac_list$gene_metadata))
filtered_enrichments = purrr::map(cluster_enrichments, ~dplyr::left_join(.,unique_motifs, by = "motif_id") %>%
                                   dplyr::semi_join(expressed_motifs, by = "motif_id") %>%
                                   dplyr::arrange(-OR_log2) )
motif_enrichment_df = ldply(filtered_enrichments, .id = "cluster_name") %>% tbl_df()

#Make a heatmap of motif enrichments
#interesting_tfs = c("FOS","BATF3","POU2F1","NFKB1","IRF1","STAT1", "MAFB", "MEF2A")
interesting_tfs = c("FOS","NFKB1","IRF1","STAT1", "MAFB", "MEF2A", "SPI1")
tf_name_casual = data_frame(new_name = c("AP-1","NF-kB","IRF","STAT1", "MAFB","MEF2", "PU.1"), tf_name = interesting_tfs) %>%
  dplyr::mutate(new_name = factor(new_name, levels = rev(new_name)))

#Filter enrichments by TF
selected_enrichments = dplyr::filter(motif_enrichment_df, tf_name %in% interesting_tfs) %>%
  dplyr::left_join(tf_name_casual, by = "tf_name")

#Make a line plot of enrichments
motif_enrichment_plot = ggplot(selected_enrichments, aes(y = new_name, x = OR_log2, xmin = ci_lower_log2, xmax = ci_higher_log2)) + 
  facet_grid(cluster_name~.) + 
  geom_point() + 
  geom_errorbarh(aes(height = 0)) +
  xlab("Log2 fold enrichment") +
  ylab("TF motif name") + 
  theme_light() +
  scale_x_continuous(expand = c(0, 0), limits = c(-4,4)) +
  theme(legend.key = element_blank()) + 
  theme(panel.margin = unit(0.2, "lines")) + 
  geom_vline(aes(xintercept = 0), size = 0.3)

ggsave("figures/main_figures/ATAC_motif_enrichment_in_clusters.pdf", plot = motif_enrichment_plot, width = 4, height = 6)



#### Perform motif enrichemnt against a background of promotoer sequences #####
fimo_promoter_hits = readr::read_delim("results/ATAC/FIMO_CISBP_promoter_matches.txt.gz", delim = "\t", col_types = c("cciicddcc"), 
                              col_names = c("motif_id","seq_name","start","end","strand","score","p_value","dummy","matched_seq"), skip = 1)

#Calculate enrichments
width_matrix = dplyr::mutate(atac_list$gene_metadata, width = end - start)
enrichment = fimoFisherTest(fimo_promoter_hits, fimo_hits_clean, bg_seq_length = 21350*2000, fg_seq_length = sum(width_matrix$width))

#Filter by gene expression
enrichment = enrichment %>%
  dplyr::left_join(unique_motifs, by = "motif_id") %>%
  dplyr::semi_join(expressed_motifs, by = "motif_id") %>%
  dplyr::arrange(-fold_enrichment)

#Import motifs from disk
cisbp_pwm_list = readRDS("results/ATAC/cisBP_PWMatrixList.rds")
cisbp_pfm_list = readRDS("results/ATAC/cisBP_PFMatrixList.rds")

#Make seqLogos for the enriched motifs
enriched_motifs = dplyr::filter(enrichment, OR > 1.2)
motif_list = as.list(cisbp_pfm_list)[enriched_motifs$motif_id]
names(motif_list) = paste(enriched_motifs$tf_name, enriched_motifs$motif_id, sep = "::")
pwm_list = purrr::map(motif_list, ~Matrix(.)) %>% purrr::map(~seqLogo::makePWM(./colSums(.)))
saveMotifListLogos(pwm_list, "results/ATAC/motif_analysis/all_enriched/")

#Select subset of motifs from each family
motif_ids = c("M2278_1.02", "M6119_1.02", "M4427_1.02","M1882_1.02", "M6308_1.02",
              "M4635_1.02", "M1925_1.02", "M2268_1.02", "M1928_1.02","M6331_1.02",
              "M5292_1.02", "M1955_1.02","M6423_1.02", "M1884_1.02","M6457_1.02",
              "M4444_1.02","M1917_1.02","M5302_1.02","M6313_1.02","M2275_1.02","M2277_1.02")
selected_motifs = dplyr::filter(enriched_motifs, motif_id %in% motif_ids) %>% 
  dplyr::mutate(p_nominal = ifelse(fisher_pvalue == 0, 1e-300, fisher_pvalue)) %>%
  dplyr::mutate(tf_name = factor(tf_name, levels = rev(tf_name)))
write.table(selected_motifs, "results/ATAC/motif_analysis/cisBP_selected_enriched_motifs.txt", 
            row.names = FALSE, quote = FALSE)

#Make a bg enrichment plot
background_enrichment_plot = ggplot(selected_motifs, aes(y = tf_name, x = OR_log2, xmin = ci_lower_log2, xmax = ci_higher_log2)) + 
  geom_point() + 
  geom_errorbarh(aes(height = 0)) +
  xlab("Log2 fold enrichment") +
  ylab("TF motif name") + 
  theme_light() +
  scale_x_continuous(expand = c(0, 0), limits = c(0,3)) +
  theme(legend.key = element_blank()) + 
  theme(panel.margin = unit(0.2, "lines"))
ggsave("figures/supplementary/ATAC_motif_enrichment_in_all_peaks.pdf", plot = background_enrichment_plot, width = 4, height = 6)


#Save selected motifs to disk
motif_list = as.list(cisbp_pfm_list)[selected_motifs$motif_id]
names(motif_list) = paste(selected_motifs$tf_name, selected_motifs$motif_id, sep = "-")
pwm_list = purrr::map(motif_list, ~Matrix(.)) %>% purrr::map(~seqLogo::makePWM(./colSums(.)))
saveMotifListLogos(pwm_list, "results/ATAC/motif_analysis/selected_motifs/", width = 8, height = 4)

#Calculate similarity matrix of PWMs
motif_pwm_list = cisbp_pwm_list[enriched_motifs$motif_id]
names(motif_pwm_list) = enriched_motifs$tf_name
sim_matrix = PWMSimilarityMatrix(motif_pwm_list, method = "Pearson")
heatmap.2(sim_matrix)



