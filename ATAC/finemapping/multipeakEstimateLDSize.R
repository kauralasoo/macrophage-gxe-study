library(ggplot2)

#Import credible sets from disk and add overlapping peaks
credible_sets_df = importCredibleSets("results/ATAC/QTLs/rasqual_credible_sets.rds", atac_list$gene_metadata)

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = purrr::map(min_pvalues_list, ~dplyr::filter(., p_eigen < fdr_thresh))

#Import Master-dependent pairs
result_list = readRDS("results/ATAC/QTLs/qtl_peak_type_assignment.rds")
master_dependent_pairs = result_list$dependents$unique_masters

#Identify master QTLs
naive_masters = dplyr::semi_join(min_pvalues_hits$naive, master_dependent_pairs, by = c("gene_id" = "master_id")) %>%
  dplyr::transmute(gene_id, is_master = TRUE)
  
#Estimate LD blocks for each caQTL
ld_block_df = credible_sets_df$naive %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::summarize(snp_count = length(snp_id), ld_region = max(pos) - min(pos)) %>%
  dplyr::left_join(naive_masters, by = "gene_id") %>%
  dplyr::mutate(is_master = ifelse(is.na(is_master), FALSE, is_master))



snp_count = ggplot(dplyr::filter(ld_block_df, snp_count > 0), aes(x = snp_count, color = is_master)) + 
  geom_density() +
  theme_light()
ggsave("figures/supplementary/caQTl_master_snp_count.pdf", plot = snp_count, width = 4, height = 3)
region_length = ggplot(dplyr::filter(ld_block_df, snp_count > 0), aes(x = ld_region, color = is_master)) + 
  geom_density() +
  theme_light()
ggsave("figures/supplementary/caQTl_master_region_length.pdf", plot = region_length, width = 4, height = 3)
