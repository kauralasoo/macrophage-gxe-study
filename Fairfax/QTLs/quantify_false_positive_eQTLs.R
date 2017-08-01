library("dplyr")
library("purrr")
library("ggplot2")

#Import minimal QTL pvalues
qtl_min_pvalues = readRDS("results/Fairfax/fairfax_qtl_min_pvalues.rds")

#Indentify hits with each sample size
full_hits = purrr::map(qtl_min_pvalues$full, ~dplyr::filter(., p_fdr < 0.1))
shared_hits = purrr::map(qtl_min_pvalues$shared, ~dplyr::filter(., p_fdr < 0.1))
shared84_hits = purrr::map(qtl_min_pvalues$shared_84, ~dplyr::filter(., p_fdr < 0.1))
shared42_hits = purrr::map(qtl_min_pvalues$shared_42, ~dplyr::filter(., p_fdr < 0.1))

#Calculate overlaps
shared84_overlap = purrr::map2(shared_hits, shared84_hits, ~dplyr::data_frame(overlap = length(intersect(.x$group_id, .y$group_id))/length(.y$group_id))) %>%
  purrr::map_df(identity, .id = "condition_name") %>%
  dplyr::mutate(sample_size = 84)
shared42_overlap = purrr::map2(shared_hits, shared42_hits, ~dplyr::data_frame(overlap = length(intersect(.x$group_id, .y$group_id))/length(.y$group_id))) %>%
  purrr::map_df(identity, .id = "condition_name") %>%
  dplyr::mutate(sample_size = 42)

olaps = dplyr::bind_rows(shared42_overlap, shared84_overlap)
write.table(olaps, "figures/supplementary/Fairfax_subsample_qtl_sharing.txt", sep = "\t", quote = FALSE, row.names = FALSE)

