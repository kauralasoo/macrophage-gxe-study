library("Biostrings")
library("TFBSTools")
library("dplyr")
library("devtools")
library("purrr")
library("ggplot2")
load_all("../seqUtils/")

#Import snp info and identify indels
snp_info = importVariantInformation("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
indels = dplyr::filter(snp_info, is_indel == TRUE) %>% dplyr::select(snp_id, is_indel)

#Import putative causal peaks/variants
unique_peaks = readRDS("results/ATAC/QTLs/qtl_peak_type_assignment.rds")$unique_masters
unique_peaks_filtered = dplyr::filter(unique_peaks$lead_snps, overlap_snp_count <= 3)
unique_overlapping_snps = dplyr::semi_join(unique_peaks$lead_credible_sets, unique_peaks_filtered, by = "gene_id") %>% 
  dplyr::filter(gene_id == overlap_peak_id) %>%
  dplyr::select(gene_id, snp_id)  %>%
  dplyr::left_join(indels, by = "snp_id") %>%
  dplyr::mutate(is_indel = ifelse(is.na(is_indel),FALSE, TRUE))

#Exluced peaks with indels from the analysis
indel_peaks = unique_overlapping_snps %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::summarise(has_indel = max(is_indel)) %>% 
  dplyr::filter(has_indel == 1)
unique_peaks_no_indels = dplyr::anti_join(unique_overlapping_snps, indel_peaks, by = "gene_id")

#Import variable chromatin QTLs
variable_qtls = readRDS("results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")
variable_peaks = rbind(dplyr::select(variable_qtls$appear, gene_id) %>% unique(), dplyr::select(variable_qtls$disappear, gene_id) %>% unique()) %>%
  dplyr::semi_join(unique_peaks_no_indels, by = "gene_id")

#Idenify peaks with with no interactions
interaction_df = readRDS("results/ATAC/QTLs/rasqual_interaction_results.rds")
no_interaction = dplyr::filter(interaction_df, p_nominal > 0.5)

#Import motif matches
motif_colnames = c("gene_id","snp_id","snp_count",".row","strand","motif_id","ref_abs_score","ref_rel_score",
                   "ref_match","alt_abs_score","alt_rel_score","alt_match","rel_diff","max_rel_score")
motif_disruptions = readr::read_delim("results/ATAC/motif_analysis/motif_disruption.txt.gz", delim = "\t", col_names = motif_colnames) %>%
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
appear_cluster_members = dplyr::group_by(variable_qtls$appear, new_cluster_id) %>% 
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

#Combine QTL clusters into larger groups to increase power
cluster_groups = data_frame(new_cluster_id = c(1:6), 
                            cluster_group_id = c("1","2-3","2-3","4","5-6","5-6"))
grouped_clusters = dplyr::left_join(appear_disruptions, cluster_groups, by = "new_cluster_id") %>%
  group_by(cluster_group_id, tf_name) %>% 
  dplyr::summarise(cluster_size = sum(cluster_size), cluster_disruption = sum(cluster_disruption)) %>%
  dplyr::rename(new_cluster_id = cluster_group_id) %>%
  ungroup()

#Calculate fold enrichment
relative_enrichment = dplyr::left_join(baseline_enrichment, grouped_clusters, by = "tf_name") %>% 
  dplyr::arrange(new_cluster_id) %>%
  purrr::by_row(., ~subsampleFisherTest(.$cluster_disruption, .$cluster_size, 
              .$baseline_disruption, .$baseline_peak_count), .collate = "rows") %>%
  dplyr::mutate(fraction_peaks = round(cluster_disruption/cluster_size,3)) %>%
  dplyr::mutate(fraction_baseline = round(baseline_disruption/baseline_peak_count,3)) %>%
  arrange(new_cluster_id, -OR_log2) %>%
  dplyr::left_join(motif_names, by = "tf_name") %>%
  #dplyr::filter(tf_name %in% c("IRF1","IRF8","RELA","NFKB1","FOS","STAT1","SPI1", "CEBPB")) %>%
  dplyr::filter(tf_name %in% c("IRF1","RELA","FOS","STAT1","SPI1")) %>%
  dplyr::mutate(tf_name = factor(tf_name, levels = rev(c("IRF1","IRF8","RELA","NFKB1","FOS","STAT1","SPI1", "CEBPB")))) %>%
  censorLog2OR()

#Add new motif names
motif_map = data_frame(tf_name = c("IRF1","RELA","FOS","STAT1","SPI1"),
                       motif_name = c("IRF","NF-kB","AP-1","STAT1","PU.1")) %>%
  dplyr::mutate(motif_name = factor(motif_name, levels = rev(motif_name)))
relative_enrichment_renamed = dplyr::left_join(relative_enrichment, motif_map, by = "tf_name")
saveRDS(relative_enrichment_renamed, "results/ATAC/motif_analysis/caQTL_clusters_enriced_motifs.rds")

#Make a lineplot of enrichments
plot = ggplot(relative_enrichment_renamed, aes(y = motif_name, x = OR_log2, xmin = ci_lower_log2, xmax = ci_higher_log2)) + 
  facet_grid(new_cluster_id~.) + 
  geom_point() + 
  geom_errorbarh(aes(height = 0)) +
  xlab("Log2 fold enrichment") +
  ylab("TF motif name") + 
  theme_light() +
  scale_x_continuous(expand = c(0, 0), limits = c(-4,4)) +
  theme(legend.key = element_blank()) + 
  theme(panel.margin = unit(0.2, "lines"))  + 
  geom_vline(aes(xintercept = 0), size = 0.3)
ggsave("figures/supplementary/caQTL_all_clusters_disrupted_motifs.pdf", plot = plot, height = 5.5, width = 3)


#Make a lineplot for Salmonella and IFNg specific clusters
filtered_df = dplyr::filter(relative_enrichment_renamed, new_cluster_id %in% c("2-3", "5-6")) %>% 
  dplyr::mutate(cluster_name = ifelse(new_cluster_id == "2-3", "Salmonella-specific", "IFNg-specific")) %>%
  dplyr::mutate(cluster_name = factor(cluster_name, levels = c("Salmonella-specific", "IFNg-specific")))

compact_enrichment_plot = ggplot(filtered_df, aes(y = motif_name, x = OR_log2, xmin = ci_lower_log2, xmax = ci_higher_log2)) + 
  facet_wrap(~cluster_name, ncol = 1) + 
  geom_point() + 
  geom_errorbarh(aes(height = 0)) +
  xlab(expression(paste(Log[2], " fold enrichment"))) +
  ylab("TF motif name") + 
  theme_light() +
  scale_x_continuous(expand = c(0, 0), limits = c(-4,4)) +
  theme(legend.key = element_blank()) + 
  theme(panel.spacing = unit(0.2, "lines"))  + 
  geom_vline(aes(xintercept = 0), size = 0.3)
ggsave("figures/main_figures/caQTL_clusters_disrupted_motifs.pdf", plot = compact_enrichment_plot, height = 3.5, width = 3)



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


