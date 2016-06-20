library("devtools")
library("purrr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Import credible sets from disk
credible_sets = readRDS("results/ATAC/QTLs/rasqual_credible_sets.rds")

#Convert into data frames and add overlapping peak annotations
credible_sets_df = purrr::map(credible_sets, ~purrr::map_df(., ~dplyr::mutate(.,chr = as.character(chr))) %>%
                                dplyr::filter(chr != "X") %>%
                                addOverlappingPeaks(atac_list$gene_metadata))

#Find lead SNPs and credible sets for unique peaks
credible_sets_by_condition = map_df(credible_sets_df, identity, .id = "condition_name")


#Extract peaks that share credible sets
shared_credible_sets = purrr::map(credible_sets_df, ~dplyr::filter(.,!is.na(overlap_peak_id))) %>% 
  purrr::map_df(., identity) %>% 
  dplyr::select(gene_id, overlap_peak_id) %>% 
  unique() %>% 
  dplyr::rename(master_peak_id = gene_id)

#Extract potential master peaks
potential_master_peaks = purrr::map(credible_sets_df, ~dplyr::filter(.,gene_id == overlap_peak_id)) %>% 
  purrr::map_df(., identity)
master_peak_list = dplyr::select(potential_master_peaks, gene_id) %>% unique()

#Identify cases where there are more than one putative master peak
shared_masters = dplyr::filter(shared_credible_sets, overlap_peak_id %in% master_peak_list$gene_id) %>%
  dplyr::filter(master_peak_id %in% master_peak_list$gene_id) %>%
  dplyr::filter(master_peak_id != overlap_peak_id)
shared_master_peaks = c(shared_masters$overlap_peak_id, shared_masters$master_peak_id) %>% unique()


#### UNIQUE MASTER PEAKS ####
#Identify unique master peaks
unique_masters = dplyr::filter(master_peak_list, !(gene_id %in% shared_master_peaks))
unique_master_pairs = dplyr::semi_join(potential_master_peaks, unique_masters, by = "gene_id") %>% dplyr::select(gene_id, snp_id) %>% unique()
snp_count = dplyr::group_by(unique_master_pairs, gene_id) %>% dplyr::summarise(snp_count = length(snp_id))
unique_masters_counted = dplyr::left_join(unique_master_pairs, snp_count, by = "gene_id")

#Find lead SNPs for unique masters accross conditions
snp_stats = dplyr::left_join(unique_masters_counted, credible_sets_by_condition, by = c("gene_id","snp_id"))
unique_lead_snps = group_by(snp_stats, gene_id) %>% 
  dplyr::arrange(p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::select(gene_id, snp_id, snp_count, condition_name, R2, chr, pos) %>%
  dplyr::ungroup()

#Find out some additional potential multi-peak QTLs
r2_filtered = filterHitsR2(dplyr::mutate(unique_lead_snps, gene_id = "gene"), vcf_file$genotypes)
overlapping_snps = dplyr::anti_join(unique_lead_snps, r2_filtered, by = "snp_id")

#For each SNP, find the other overlapping SNPs
genotype_matrix = t(vcf_file$genotypes[unique_lead_snps$snp_id,])
r2 = cor(genotype_matrix, use = "pairwise.complete.obs")^2
olap_snps = r2[overlapping_snps$snp_id,]

results = list()
for(snp_id in overlapping_snps$snp_id) {
  olap_ids = names(which(olap_snps[snp_id, ] > 0.8))
  df = data_frame(master_snp_id = snp_id, overlap_snp_id = olap_ids)
  results[[snp_id]] = df
}
shared_snps = map_df(results, identity) %>% dplyr::filter(master_snp_id != overlap_snp_id)
shared_snp_ids = c(shared_snps$master_snp_id, shared_snps$overlap_snp_id) %>% unique()

#Remove additional LD friends
unique_lead_snps_filtered = dplyr::filter(unique_lead_snps, !(snp_id %in% shared_snp_ids))
unique_credible_sets = dplyr::semi_join(credible_sets_by_condition, unique_lead_snps_filtered, by = c("gene_id", "condition_name"))
unique_masters_counted = dplyr::filter(unique_masters_counted, gene_id %in% unique_credible_sets$gene_id)

#Combine unique results
unique_peaks = list(peak_snp_pairs = unique_masters_counted, lead_snps = unique_lead_snps_filtered, lead_credible_sets = unique_credible_sets)
saveRDS(unique_peaks, "results/ATAC/QTLs/unique_qtl_peaks.rds")


#### Multiple master peaks per QTL ####

#Deal with lead SNPs that were not unique
peak_snp_map = dplyr::select(unique_lead_snps, gene_id, snp_id)
peaks_with_shared_lead_snps = dplyr::left_join(shared_snps, peak_snp_map, by = c("master_snp_id" = "snp_id")) %>% 
  dplyr::rename(master_peak_id = gene_id) %>% 
  dplyr::left_join(peak_snp_map, by =c("overlap_snp_id" = "snp_id")) %>% 
  dplyr::rename(overlap_peak_id = gene_id) %>% 
  dplyr::select(master_peak_id, overlap_peak_id)
shared_masters = rbind(shared_masters, peaks_with_shared_lead_snps) %>% unique()
shared_lead_peaks = unique(c(peaks_with_shared_lead_snps$master_peak_id, peaks_with_shared_lead_snps$overlap_peak_id))

#Convert shared peaks into clusters
atac_clusters = constructClustersFromGenePairs(shared_masters, cluster_name_prefix = "ATAC_cluster_")

#Count the number of SNPs in each peak and cluster
atac_clusters_counted = dplyr::left_join(atac_clusters, potential_master_peaks, by = "gene_id") %>% 
  dplyr::group_by(gene_id, cluster_id) %>% 
  dplyr::summarise(peak_snp_count = length(snp_id)) %>% 
  dplyr::group_by(cluster_id) %>% dplyr::mutate(cluster_snp_count = sum(peak_snp_count), peak_count = length(gene_id)) %>% 
  dplyr::ungroup() %>%
  arrange(cluster_snp_count)

#Find lead SNPs for each cluster of peaks
cluster_lead_snps = dplyr::left_join(atac_clusters_counted, potential_master_peaks, by = "gene_id") %>% 
  dplyr::group_by(cluster_id) %>% 
  dplyr::arrange(p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::select(cluster_id, gene_id, snp_id, R2, peak_count, cluster_snp_count, chr, pos) %>%
  dplyr::ungroup()

#Save cluster results to disk
cluster_peaks = list(cluster_memberships = atac_clusters_counted, lead_snps = cluster_lead_snps)
saveRDS(cluster_peaks, "results/ATAC/QTLs/clustered_qtl_peaks.rds")

### Dependent peaks with significant master QTLs ####
#Find all caQTLs whose credible set overlaps some other peak, but not the peak itself
no_master_qtls = dplyr::filter(shared_credible_sets, master_peak_id != overlap_peak_id) %>% 
  dplyr::anti_join(master_peak_list, by = c("master_peak_id" = "gene_id"))

#Find all of the cases where the other peak is a unique master caQTL peak
dependent_uniq_masters = dplyr::semi_join(no_master_qtls, unique_lead_snps_filtered,by = c("overlap_peak_id" = "gene_id")) %>%
  dplyr::left_join(unique_lead_snps_filtered, by = c("overlap_peak_id" = "gene_id")) %>% 
  dplyr::rename(dependent_id = master_peak_id, master_id = overlap_peak_id)

#Find all of the cases where the other peak is belongs to a cluster of master peaks
dependent_cluster_masters = dplyr::semi_join(no_master_qtls, cluster_lead_snps,by = c("overlap_peak_id" = "gene_id")) %>%
  dplyr::left_join(cluster_lead_snps, by = c("overlap_peak_id" = "gene_id")) %>% 
  dplyr::rename(dependent_id = master_peak_id, master_id = overlap_peak_id)
res = list(unique_masters = dependent_uniq_masters, cluster_master = dependent_cluster_masters)
saveRDS(res, "results/ATAC/QTLs/depdendent_peaks.rds")

##### Do summary stats ####
total_peak_count = map_df(credible_sets_df, identity) %>% dplyr::select(gene_id) %>% unique() %>% nrow()
overlap_peak_count = purrr::map(credible_sets_df, ~dplyr::filter(.,!is.na(overlap_peak_id))) %>% 
  purrr::map_df(., identity) %>% 
  dplyr::select(gene_id) %>% unique() %>% nrow()
overlap_same_peak = dplyr::select(potential_master_peaks, gene_id) %>% unique() %>% nrow()
unique_master_peaks = length(unique_masters_counted$gene_id %>% unique)
shared_master_peaks = nrow(atac_clusters_counted)
other_is_qtl = length(unique(master_cluster$master_peak_id)) + length(unique(master_unique$master_peak_id))

#Compile stats
summary_stats = data_frame("total" = total_peak_count, "overlap_any_peak" = overlap_peak_count, "overlap_same_peak" = overlap_same_peak, 
                           "shared_masters" = shared_master_peaks, "unique_masters" = unique_master_peaks, "other_is_qtl" = other_is_qtl) %>%
  dplyr::mutate(overlap_other_peak = overlap_any_peak-overlap_same_peak, overlap_no_peak = total - overlap_any_peak) %>%
  dplyr::mutate(other_not_qtl = overlap_other_peak - other_is_qtl)
write.table(summary_stats, "results/ATAC/QTLs/properties/credible_set_peak_overlaps.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Make a plot of counts
summary_df = dplyr::select(summary_stats, unique_masters, shared_masters, other_is_qtl, other_not_qtl, overlap_no_peak) %>% 
  tidyr::gather("type", "count")
plot = ggplot(summary_df, aes(x = factor(1), fill = type, y = count)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette="Set1") +
  theme(legend.position="none")
ggsave("results/ATAC/QTLs/properties/caQTL_classification.pdf", plot = plot, width = 1.4, height = 4)


#Plot the number of peaks per cluster
cluster_size_df = dplyr::transmute(atac_clusters_counted, cluster_id, peak_count = factor(peak_count)) %>% unique()
count_plot = ggplot(cluster_size_df, aes(x = peak_count)) + geom_bar() + theme_hc() +
  xlab("Number of peaks")
ggsave("results/ATAC/QTLs/properties/number_of_peaks_per_cluster.pdf", plot = count_plot, width = 4, height = 4)

#Plot the number of SNPs per unique QTL peak
snp_count_df = dplyr::select(unique_masters_counted, gene_id, snp_count) %>% unique() %>%
  dplyr::filter(snp_count <= 10) %>%
  dplyr::mutate(snp_count = factor(snp_count))
snp_count_plot = ggplot(snp_count_df, aes(x = snp_count)) + geom_bar() + theme_hc() +
  xlab("Number of variants")
ggsave("results/ATAC/QTLs/properties/number_of_variants_per_unqiue_master.pdf", plot = snp_count_plot, width = 4, height = 4)

#Count the number of dependent peaks per master
master_dependent_pairs = rbind(dplyr::select(dependent_uniq_masters, dependent_id, master_id), dplyr::select(dependent_cluster_masters, dependent_id, master_id))
dependent_peak_count = master_dependent_pairs %>% 
  dplyr::group_by(master_id) %>% 
  dplyr::summarise(dependent_peak_count = length(dependent_id)) %>% 
  arrange(-dependent_peak_count)
dependent_plot = ggplot(dependent_peak_count, aes(x = dependent_peak_count)) + geom_bar() + theme_hc() +
  xlab("Number of dependent peaks")
ggsave("results/ATAC/QTLs/properties/number_of_dependent_peaks_per_unqiue_master.pdf", plot = dependent_plot, width = 4, height = 4)

#Plot distances between master and dependent peaks
dependent_distances = calculatePeakDistance(master_dependent_pairs, atac_data$gene_metadata)
dependent_distance_plot = ggplot(dependent_distances, aes(x = distance/1000)) + geom_histogram(binwidth = 2) + theme_hc() +
  xlab("Distance from master peak")
ggsave("results/ATAC/QTLs/properties/master_dependent_peak_distance.pdf", plot = dependent_distance_plot, width = 4, height = 4)

#Plot the number of variants per cluster
cluster_variants = ggplot(cluster_lead_snps, aes(x = cluster_snp_count)) + 
  geom_bar() + 
  scale_x_continuous(limits = c(0,50)) + 
  theme_hc() + 
  xlab("Number of variants per cluster")
ggsave("results/ATAC/QTLs/properties/number_of_variants_per_cluster.pdf", plot = cluster_variants, width = 4, height = 4)

#Calculate distances within clusters
cluster_dist = dplyr::rename(shared_masters, master_id = master_peak_id, dependent_id = overlap_peak_id) %>% 
  calculatePeakDistance(atac_list$gene_metadata)
cluster_dist_plot = ggplot(cluster_dist, aes(x = distance/1000)) + geom_histogram(binwidth = 2) + theme_hc() +
  xlab("Distance between two peaks in cluster") + 
  scale_x_continuous(limits = c(-50, 50))
ggsave("results/ATAC/QTLs/properties/within_cluster_distances.pdf", plot = cluster_dist_plot, width = 4, height = 4)


