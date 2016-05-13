library("Biostrings")
library("TFBSTools")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("purrr")
library("ggplot2")

#Import snp info and identify indels
snp_info = readr::read_delim("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz", 
                             delim = "\t", col_types = "cdccc", col_names = c("chr","pos","snp_id","ref","alt"))
snp_info = dplyr::mutate(snp_info, indel_length = pmax(nchar(alt), nchar(ref))) %>%
  dplyr::mutate(is_indel = ifelse(indel_length > 1, TRUE, FALSE))
indels = dplyr::filter(snp_info, is_indel == TRUE) %>% dplyr::select(snp_id, is_indel)

#Import putative causal peaks/variants
unique_peaks = readRDS("results/ATAC/QTLs/unique_qtl_peaks.rds")
unique_peaks_filtered = dplyr::filter(unique_peaks$peak_snp_pairs, snp_count <= 3) %>%
  dplyr::left_join(indels, by = "snp_id") %>%
  dplyr::mutate(is_indel = ifelse(is.na(is_indel),FALSE, TRUE))
#Exluced peaks with indels from the analysis
indel_peaks = unique_peaks_filtered %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::summarise(has_indel = max(is_indel)) %>% 
  dplyr::filter(has_indel == 1)
unique_peaks_no_indels = dplyr::anti_join(unique_peaks_filtered, indel_peaks, by = "gene_id")

#Import variable chromatin QTLs
variable_qtls = readRDS("results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")
variable_peaks = rbind(dplyr::select(variable_qtls$appear, gene_id) %>% unique(), dplyr::select(variable_qtls$disappear, gene_id) %>% unique()) %>%
  dplyr::semi_join(unique_peaks_no_indels, by = "gene_id")

#Idenify peaks with with no interactions
interaction_df = readRDS("results/ATAC/QTLs/rasqual_interaction_results.rds")
no_interaction = dplyr::filter(interaction_df, p_nominal > 0.5)

#Import motif matches
motif_colnames = c("gene_id","snp_id","snp_count",".row","start","strand","motif_id","ref_abs_score","ref_rel_score",
                   "ref_match","alt_abs_score","alt_rel_score","alt_match","rel_diff","max_rel_score")
motif_disruptions = readr::read_delim("results/ATAC/motif_analysis/motif_disruption.txt", delim = "\t", col_names = motif_colnames) %>%
  dplyr::anti_join(indels, by = "snp_id")

#Load TF names
mf_enriched_motifs = read.table("results/ATAC/motif_analysis/cisBP_selected_enriched_motifs.txt", header = TRUE, stringsAsFactors = FALSE)
motif_names = dplyr::select(mf_enriched_motifs, motif_id, tf_name)
motif_class = read.table("results/ATAC/motif_analysis/tf_name_class_map.txt", header = TRUE, stringsAsFactors = FALSE)
motif_names = dplyr::left_join(motif_names, motif_class, by = "tf_name")

#Extract motif dsiruption events
disruption_events = dplyr::filter(motif_disruptions, max_rel_score > 0.85, rel_diff > 0.03) %>%
  dplyr::left_join(motif_names, by = "motif_id")

#Identfy motif disruptions in all caQTLs
baseline_enrichment = dplyr::semi_join(disruption_events, unique_peaks_no_indels, by = "gene_id") %>% 
  dplyr::group_by(gene_id, tf_name) %>% 
  dplyr::arrange(-abs(rel_diff)) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::group_by(tf_name) %>% 
  dplyr::summarise(baseline_disruption = length(gene_id)) %>% 
  dplyr::arrange(-baseline_disruption) %>%
  dplyr::mutate(baseline_peak_count = length(unique(unique_peaks_no_indels$gene_id)))
write.table(baseline_enrichment, "results/ATAC/motif_analysis/baseline_disrupted_motifs.txt", sep ="\t", row.names = FALSE, quote = FALSE)

#Calculate the proportion of peaks affected
proportion_disrupted_by_motif = dplyr::mutate(baseline_enrichment, proportion = baseline_disruption/baseline_peak_count)
total_disruptions = dplyr::semi_join(disruption_events, unique_peaks_no_indels, by = "gene_id") %>% 
  dplyr::group_by(gene_id, tf_name) %>% 
  dplyr::arrange(-abs(rel_diff)) %>% 
  dplyr::filter(row_number() == 1) %>% 
  group_by(gene_id) %>% 
  summarise(motif_disruption_count = length(tf_name))
fraction_of_peaks_disrupted = nrow(total_disruptions)/ length(unique(unique_peaks_no_indels$gene_id))

#Analyse QTLs that appear
appear_cluster_members = dplyr::group_by(variable_qtls$appear, cluster_id ) %>% 
  purrr::by_slice(~dplyr::select(.,gene_id) %>% 
                    unique() %>%
                    dplyr::semi_join(unique_peaks_no_indels, by = "gene_id"), .to = "clusters")
appear_cluster_members = dplyr::mutate(appear_cluster_members, cluster_size = map(appear_cluster_members$clusters, nrow) %>% unlist())

#Count motif disruptions for Appearing QTL clusters
appear_disruptions = appear_cluster_members %>% 
  purrr::by_row(~dplyr::semi_join(disruption_events, .$clusters[[1]]) %>%
                                             dplyr::group_by(gene_id, tf_name) %>% 
                                             dplyr::arrange(-abs(rel_diff)) %>% dplyr::filter(row_number() == 1) %>%
                                             dplyr::group_by(tf_name) %>% 
                                             dplyr::summarise(cluster_disruption = length(gene_id)) %>%
                                             dplyr::left_join(dplyr::select(motif_names, tf_name), ., by = "tf_name") %>% 
                                             dplyr::mutate(cluster_disruption = ifelse(is.na(cluster_disruption), 0, cluster_disruption)) %>%
                                             dplyr::arrange(-cluster_disruption), .collate = "rows") %>%
  dplyr::select(-clusters, -.row)

#Calculate fold enrichment
relative_enrichment = dplyr::left_join(baseline_enrichment, appear_disruptions, by = "tf_name") %>% 
  dplyr::arrange(cluster_id) %>% 
  dplyr::mutate(fold_enrichment = (cluster_disruption/cluster_size)/(baseline_disruption/baseline_peak_count)) %>% 
  dplyr::mutate(fold_enrichment = round(ifelse(fold_enrichment == 0, 0.1, fold_enrichment),3)) %>%
  dplyr::mutate(l2_fold = log(fold_enrichment, 2)) %>%
  dplyr::mutate(fraction_peaks = round(cluster_disruption/cluster_size,3)) %>%
  dplyr::mutate(fraction_baseline = round(baseline_disruption/baseline_peak_count,3)) %>%
  arrange(cluster_id, -fold_enrichment) %>%
  dplyr::left_join(motif_names, by = "tf_name") %>%
  dplyr::filter(tf_name %in% c("IRF1","IRF8","RELA","NFKB1","FOS","STAT1","SPI1", "CEBPB")) %>%
  dplyr::mutate(tf_name = factor(tf_name, levels = rev(c("IRF1","IRF8","RELA","NFKB1","FOS","STAT1","SPI1", "CEBPB"))))

#Calculate enrichment p-values
appear_enrihced_motifs = dplyr::mutate(relative_enrichment, p_enriched = phyper(cluster_disruption -1, baseline_disruption, baseline_peak_count - baseline_disruption, cluster_size, lower.tail = FALSE)) %>% 
  dplyr::mutate(fdr_enriched = p.adjust(p_enriched, "fdr")) %>%
  dplyr::mutate(p_depleted = phyper(cluster_disruption, baseline_disruption, baseline_peak_count - baseline_disruption, cluster_size, lower.tail = TRUE)) %>%
  dplyr::mutate(fdr_depleted = p.adjust(p_depleted, "fdr")) %>%
  dplyr::mutate(is_significant = ifelse(fdr_enriched > 0.1, FALSE, TRUE))
saveRDS(appear_enrihced_motifs, "results/ATAC/motif_analysis/caQTL_clusters_enriced_motifs.rds")

#Identify the fraction of peaks that disrupt enriched motifs
fraction_disrupting_enriched_motifs = dplyr::filter(appear_enrihced_motifs, is_significant == TRUE) %>% 
  dplyr::select(tf_name, cluster_id, cluster_size, cluster_disruption, fraction_baseline, fraction_peaks, fold_enrichment, p_enriched)
write.table(fraction_disrupting_enriched_motifs, "results/ATAC/motif_analysis/fraction_disrupting_enriched_motifs.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#Make a plot of enriched motifs
plot = ggplot(appear_enrihced_motifs, aes(y = tf_name, x = l2_fold, color = is_significant, size = -log10(p_enriched))) + 
  facet_grid(cluster_id~.) + 
  geom_point() + 
  scale_size(range = c(2,8)) + 
  xlab("Log2 fold enrichment") +
  ylab("TF motif name")
ggsave("results/ATAC/motif_analysis/caQTL_clusters_enriched_motifs.pdf", plot = plot, height = 8, width = 6)




#Count motif disruptions for disappearing caQTL clusters
disappear_cluster_members = dplyr::group_by(variable_qtls$disappear, cluster_id ) %>% 
  purrr::by_slice(~dplyr::select(.,gene_id) %>% 
                    unique() %>%
                    dplyr::semi_join(unique_peaks_no_indels, by = "gene_id"), .to = "clusters")
disappear_cluster_members = dplyr::mutate(disappear_cluster_members, cluster_size = map(disappear_cluster_members$clusters, nrow) %>% unlist())


disappear_disruptions = disappear_cluster_members %>% 
  purrr::by_row(~dplyr::semi_join(disruption_events, .$clusters[[1]]) %>%
                  dplyr::group_by(gene_id, tf_name) %>% 
                  dplyr::arrange(-abs(rel_diff)) %>% dplyr::filter(row_number() == 1) %>%
                  dplyr::group_by(tf_name) %>% 
                  dplyr::summarise(cluster_disruption = length(gene_id)) %>% 
                  dplyr::arrange(-cluster_disruption), .collate = "rows") %>%
  dplyr::select(-clusters, -.row)

relative_enrichment = dplyr::left_join(baseline_enrichment, disappear_disruptions, by = "tf_name") %>% 
  dplyr::arrange(cluster_id) %>% 
  dplyr::mutate(fold_enrichment = (cluster_disruption/cluster_size)/(baseline_disruption/baseline_peak_count)) %>% 
  arrange(cluster_id, -fold_enrichment)

disappear_enrihced_motifs = dplyr::mutate(relative_enrichment, p_nominal = phyper(cluster_disruption -1, baseline_disruption, baseline_peak_count - baseline_disruption, cluster_size, lower.tail = FALSE)) %>% 
  dplyr::mutate(p_fdr = p.adjust(p_nominal, "fdr")) %>% dplyr::filter(p_fdr < 0.1)
#RESULT: no enrichment


#Analyse all disappearing QTLs together
disappear_qtls = dplyr::filter(variable_qtls$disappear, cluster_id == 2) %>% dplyr::select(gene_id) %>% unique() %>% dplyr::semi_join(unique_peaks_no_indels, by = "gene_id")
disappear_disruptions = dplyr::semi_join(disruption_events, disappear_qtls) %>%
  dplyr::group_by(gene_id, tf_name) %>% 
  dplyr::arrange(-abs(rel_diff)) %>% dplyr::filter(row_number() == 1) %>%
  dplyr::group_by(tf_name) %>% 
  dplyr::summarise(cluster_disruption = length(gene_id)) %>% 
  dplyr::arrange(-cluster_disruption) %>% 
  dplyr::mutate(cluster_size = nrow(disappear_qtls))
dis_relative_enrichment = dplyr::left_join(baseline_enrichment, disappear_disruptions, by = "tf_name") %>% 
  dplyr::mutate(fold_enrichment = (cluster_disruption/cluster_size)/(baseline_disruption/baseline_peak_count)) %>% 
  dplyr::mutate(fold_enrichment = ifelse(fold_enrichment == 0, 0.1, fold_enrichment)) %>%
  dplyr::mutate(l2_fold = log(fold_enrichment, 2)) %>%
  arrange(-fold_enrichment) %>%
  dplyr::mutate(p_nominal = phyper(cluster_disruption -1, baseline_disruption, baseline_peak_count - baseline_disruption, cluster_size, lower.tail = FALSE)) %>% 
  dplyr::mutate(p_fdr = p.adjust(p_nominal, "fdr"))

#RESULT: no enrichment


