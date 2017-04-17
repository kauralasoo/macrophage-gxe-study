library("devtools")
library("purrr")
library("dplyr")
library("ggplot2")
library("scales")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = purrr::map(min_pvalues_list, ~dplyr::filter(., p_eigen < fdr_thresh))

#Import credible sets from disk and add overlapping peaks
credible_sets_df = importCredibleSets("results/ATAC/QTLs/rasqual_credible_sets.rds", atac_list$gene_metadata)

#Extract all peak-peak pairs that are connected by credible sets
shared_credible_sets = credibleSetPeakOverlaps(credible_sets_df)

#Identify all potential master peaks
all_master_peaks = dplyr::filter(shared_credible_sets, master_peak_id == overlap_peak_id) %>% 
  dplyr::select(master_peak_id)

#Identify cases where there are more than one putative master peak
ambiguous_master_peaks = identifyAmbiguousMasters(shared_credible_sets, all_master_peaks)

#Identify unique master peaks and their lead SNPs
unique_master_peaks = identifyUniqueMasters(all_master_peaks, ambiguous_master_peaks, credible_sets_df)

#Identify additional linked lead SNP pairs (R2 > 0.8) from the unique master peaks.
linked_lead_snps = identifyLinkedLeadSNPs(unique_master_peaks, vcf_file$genotypes)
linked_lead_snp_ids = as.vector(unlist(linked_lead_snps))

#Remove additional LD friends
unique_master_peaks_filtered = dplyr::filter(unique_master_peaks, !(snp_id %in% linked_lead_snp_ids))
unique_credible_sets = purrr::map_df(credible_sets_df, identity, .id = "condition_name") %>%
  dplyr::semi_join(unique_master_peaks_filtered, by = c("gene_id", "condition_name"))

#Combine unique results
unique_peaks = list(lead_snps = unique_master_peaks_filtered, lead_credible_sets = unique_credible_sets)


#### Identify ambiguous master caQTL clusters ####
shared_masters = dplyr::semi_join(shared_credible_sets, ambiguous_master_peaks, by = "master_peak_id") %>% 
  dplyr::semi_join(ambiguous_master_peaks, by = c("overlap_peak_id" = "master_peak_id")) %>% 
  dplyr::filter(master_peak_id != overlap_peak_id)

#Deal with lead SNPs that were not unique
peak_snp_map = dplyr::select(unique_master_peaks, gene_id, snp_id)
peaks_with_shared_lead_snps = dplyr::left_join(linked_lead_snps, peak_snp_map, by = c("master_snp_id" = "snp_id")) %>% 
  dplyr::rename(master_peak_id = gene_id) %>% 
  dplyr::left_join(peak_snp_map, by =c("overlap_snp_id" = "snp_id")) %>% 
  dplyr::rename(overlap_peak_id = gene_id) %>% 
  dplyr::select(master_peak_id, overlap_peak_id)

#Cluster ambiguous master peaks
ambiguous_masters_clustered = rbind(shared_masters, peaks_with_shared_lead_snps) %>% 
  unique() %>%
  clusterAmbiguousMasters(credible_sets_df)


### Dependent peaks with significant master QTLs ####
#Find all caQTLs whose credible set overlaps some other peak, but not the peak itself
not_master_qtls = dplyr::filter(shared_credible_sets, master_peak_id != overlap_peak_id) %>% 
  dplyr::anti_join(all_master_peaks, by = "master_peak_id")

#Find all of the cases where the other peak is a unique master caQTL peak
dependent_uniq_masters = dplyr::semi_join(not_master_qtls, unique_master_peaks_filtered, by = c("overlap_peak_id" = "gene_id")) %>%
  dplyr::left_join(unique_master_peaks_filtered, by = c("overlap_peak_id" = "gene_id")) %>% 
  dplyr::rename(dependent_id = master_peak_id, master_id = overlap_peak_id)

#Find all of the cases where the other peak is belongs to a cluster of master peaks
dependent_cluster_masters = dplyr::semi_join(not_master_qtls, ambiguous_masters_clustered$lead_snps,
                                             by = c("overlap_peak_id" = "gene_id")) %>%
  dplyr::left_join(ambiguous_masters_clustered$lead_snps, by = c("overlap_peak_id" = "gene_id")) %>% 
  dplyr::rename(dependent_id = master_peak_id, master_id = overlap_peak_id)
dependent_list = list(unique_masters = dependent_uniq_masters, cluster_master = dependent_cluster_masters)

#Compile all results together
result_list = list(unique_masters = unique_peaks, 
                   ambiguous_masters = ambiguous_masters_clustered, 
                   dependents = dependent_list)
saveRDS(result_list, "results/ATAC/QTLs/qtl_peak_type_assignment.rds")
result_list = readRDS("results/ATAC/QTLs/qtl_peak_type_assignment.rds")




##### Do summary stats ####
total_peak_count = purrr::map_df(credible_sets_df, identity) %>% dplyr::select(gene_id) %>% unique() %>% nrow()
overlap_peak_count = purrr::map(credible_sets_df, ~dplyr::filter(.,!is.na(overlap_peak_id))) %>% 
  purrr::map_df(., identity) %>% 
  dplyr::select(gene_id) %>% unique() %>% nrow()
overlap_same_peak = all_master_peaks %>% nrow()
unique_master_peaks = nrow(result_list$unique_masters$lead_snps)
shared_master_peaks = nrow(result_list$ambiguous_masters$cluster_memberships)
other_is_qtl = length(unique(result_list$dependents$unique_masters$dependent_id)) + 
  length(unique(result_list$dependents$cluster_master$dependent_id))
master_count = length(unique(result_list$dependents$unique_masters$master_id)) + 
  length(unique(result_list$dependents$cluster_master$master_id))
  

#Compile stats
summary_stats = data_frame("total" = total_peak_count, "overlap_any_peak" = overlap_peak_count, "overlap_same_peak" = overlap_same_peak, 
                           "shared_masters" = shared_master_peaks, "unique_masters" = unique_master_peaks, "other_is_qtl" = other_is_qtl) %>%
  dplyr::mutate(overlap_other_peak = overlap_any_peak-overlap_same_peak, overlap_no_peak = total - overlap_any_peak) %>%
  dplyr::mutate(other_not_qtl = overlap_other_peak - other_is_qtl)
saveRDS(summary_stats, "results/ATAC/QTLs/qtl_peak_type_summary_stats.rds")
summary_stats = readRDS("results/ATAC/QTLs/qtl_peak_type_summary_stats.rds")

##### Count the number of master caQTLs #####
summary_df = dplyr::select(summary_stats, unique_masters, shared_masters, other_is_qtl, other_not_qtl, overlap_no_peak) %>% 
  tidyr::gather("type", "count") %>%
  dplyr::mutate(fraction = count/sum(count)) %>% 
  dplyr::mutate(caqtl_type = factor(c("Same peak", "Multiple peaks", "Other QTL", "Other not QTL", "No peak"),
                                    levels = rev(c("Same peak","Other QTL", "Multiple peaks", "Other not QTL", "No peak"))))

#Make a plot of fractions
cs_plot = ggplot(summary_df, aes(x = caqtl_type, y = fraction, label = count)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  scale_y_continuous(labels=percent, breaks=pretty_breaks(n=4), limits = c(0,0.6)) +
  theme_light() +
  xlab("Credible set contains") +
  ylab("Fraction of caQTLs") + 
  coord_flip() +
  geom_text(nudge_y = 0.1) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

ggsave("figures/main_figures/caQTL_credible_set_contents.pdf", cs_plot, width = 1.7, height = 3)


##### Plot the number of SNPs per unique QTL peak ####
#Count the numbers of peaks with different numbers of SNPs in them
snp_count_df = dplyr::select(result_list$unique_masters$lead_snps, gene_id, overlap_snp_count) %>% unique() %>% 
  dplyr::mutate(overlap_snp_count = ifelse(overlap_snp_count > 9, 10, overlap_snp_count)) %>% 
  dplyr::group_by(overlap_snp_count) %>% 
  dplyr::summarise(peak_count = length(overlap_snp_count)) %>% 
  dplyr::mutate(overlap_snp_count = ifelse(overlap_snp_count == 10, ">9", as.character(overlap_snp_count))) %>%
  dplyr::mutate(overlap_snp_count = factor(overlap_snp_count, overlap_snp_count))

#Make a plot
snp_count_plot = ggplot(snp_count_df, aes(x = overlap_snp_count, y = peak_count)) + 
  geom_bar(stat = "identity") + 
  theme_light() +
  ylab("Number of master regions") +
  xlab("Number of variants in region") +
  scale_y_continuous(breaks=pretty_breaks(n=6))
ggsave("figures/main_figures/caQTL_number_of_variants_per_unqiue_master.pdf", plot = snp_count_plot, width = 3, height = 3)


##### Count the number of dependent peaks per master ####
master_dependent_pairs = rbind(dplyr::select(result_list$dependents$unique_masters, dependent_id, master_id), 
                               dplyr::select(result_list$dependents$cluster_master, dependent_id, master_id))
dependent_peak_count = master_dependent_pairs %>% 
  dplyr::group_by(master_id) %>% 
  dplyr::summarise(dependent_peak_count = length(dependent_id)) %>% 
  arrange(-dependent_peak_count)
dependent_plot = ggplot(dependent_peak_count, aes(x = dependent_peak_count)) + geom_bar() + theme_light() +
  xlab("Number of dependent regions") +
  ylab("Master region count")
ggsave("figures/main_figures/caQTL_number_of_dependent_peaks_per_unqiue_master.pdf", plot = dependent_plot, width = 3, height = 3)

#Plot distances between master and dependent peaks
dependent_distances = calculatePeakDistance(master_dependent_pairs, atac_list$gene_metadata)
dependent_distance_plot = ggplot(dependent_distances, aes(x = distance/1000)) + geom_histogram(binwidth = 2) + theme_light() +
  xlab("Distance from master region (kb)") +
  ylab("Dependent region count")
ggsave("figures/main_figures/caQTL_master_dependent_peak_distance.pdf", plot = dependent_distance_plot, width = 3, height = 3)




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
cluster_dist_plot = ggplot(cluster_dist, aes(x = distance/1000)) + geom_histogram(binwidth = 2) + theme_light() +
  xlab("Distance between two peaks in cluster") + 
  scale_x_continuous(limits = c(-50, 50))
ggsave("results/ATAC/QTLs/properties/within_cluster_distances.pdf", plot = cluster_dist_plot, width = 4, height = 4)


