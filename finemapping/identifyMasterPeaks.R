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

#Put together all credible sets
map_df(credible_sets_df, identity) %>% dplyr::select(gene_id, snp_id) %>% unique()


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

#Identify unique master peaks
unique_masters = dplyr::filter(master_peak_list, !(gene_id %in% shared_master_peaks))
unique_master_pairs = dplyr::semi_join(potential_master_peaks, unique_masters, by = "gene_id") %>% dplyr::select(gene_id, snp_id) %>% unique()
snp_count = dplyr::group_by(unique_master_pairs, gene_id) %>% dplyr::summarise(snp_count = length(snp_id))
unique_masters_counted = dplyr::left_join(unique_master_pairs, snp_count, by = "gene_id")

#Do summary stats
total_peak_count = map_df(credible_sets_df, identity) %>% dplyr::select(gene_id) %>% unique() %>% nrow()
overlap_peak_count = purrr::map(credible_sets_df, ~dplyr::filter(.,!is.na(overlap_peak_id))) %>% 
  purrr::map_df(., identity) %>% 
  dplyr::select(gene_id) %>% unique() %>% nrow()
overlap_same_peak = dplyr::select(potential_master_peaks, gene_id) %>% unique() %>% nrow()
unique_master_peaks = length(unique_masters_counted$gene_id %>% unique)

total_peak_count
overlap_peak_count
overlap_same_peak
unique_master_peaks

#Find lead SNPs and credible sets for unique peaks
credible_sets_by_condition = map_df(credible_sets_df, identity, .id = "condition_name")

#Lead SNPs
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
unique_lead_snps = dplyr::filter(unique_lead_snps, !(snp_id %in% shared_snp_ids))
unique_masters_counted = dplyr::filter(unique_masters_counted, gene_id %in% unique_credible_sets$gene_id)

#Credible sets
unique_credible_sets = dplyr::semi_join(credible_sets_by_condition, unique_lead_snps, by = c("gene_id", "condition_name"))

#Combine unique results
unique_peaks = list(peak_snp_pairs = unique_masters_counted, lead_snps = unique_lead_snps, lead_credible_sets = unique_credible_sets)
saveRDS(unique_peaks, "results/ATAC/QTLs/unique_qtl_peaks.rds")



#### Multiple master peaks per QTL ####

#Identify cases where there are more than one putative master peak
shared_masters = dplyr::filter(shared_credible_sets, overlap_peak_id %in% master_peak_list$gene_id) %>%
  dplyr::filter(master_peak_id %in% master_peak_list$gene_id) %>%
  dplyr::filter(master_peak_id != overlap_peak_id)
shared_master_peaks = c(shared_masters$overlap_peak_id, shared_masters$master_peak_id) %>% unique()

#Convert shared peaks into clusters
graph = igraph::graph_from_data_frame(shared_masters, directed = FALSE)
clusters = igraph::clusters(graph)
atac_clusters = data_frame(gene_id = names(clusters$membership), cluster_number = clusters$membership) %>% 
  arrange(cluster_number) %>% 
  dplyr::mutate(cluster_id = paste0("ATAC_cluster_",cluster_number)) %>% 
  dplyr::select(gene_id, cluster_id)
peak_counts = dplyr::group_by(atac_clusters, cluster_id) %>% dplyr::summarise(peak_count = length(gene_id))
atac_clusters = dplyr::left_join(atac_clusters, peak_counts, by = "cluster_id")

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


