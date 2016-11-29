library("ggplot2")

gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt",
                                col_names = c("trait","file_name"))

#Name all files
file_names = as.list(paste0("results/SL1344/coloc/GARFIELD/garfield.perm.",gwas_stats_labeled$file_name, ".out"))
names(file_names) = gwas_stats_labeled$trait

#Import enrichments
enrich_df = purrr::map_df(file_names, ~readr::read_delim(., delim = " "), .id = "trait")
filtered_enrichments = dplyr::filter(enrich_df, PThresh == 1e-5) %>%
  dplyr::filter(!(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "RA_2012", "CD_2012", "T1D_2", "T2D_1")))

#Calculate max fold enrichments
mean_fe_df = dplyr::group_by(filtered_enrichments, trait) %>% 
  dplyr::mutate(mean_FE = mean(FE)) %>% arrange(mean_FE) %>%
  dplyr::ungroup()
ranked_enrichments = dplyr::mutate(filtered_enrichments, trait = factor(trait, levels = unique(mean_fe_df$trait)))

#Make plots
enrichment_plot = ggplot(ranked_enrichments, aes(x = log(FE,2), y = trait, color = Annotation)) + 
  geom_point() +
  scale_x_continuous(limits = c(0, 6)) + 
  theme_light()
ggsave("figures/main_figures/grafield_enrichment.pdf", plot = enrichment_plot, width = 6, height = 7)

