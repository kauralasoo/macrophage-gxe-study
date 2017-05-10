library("ggplot2")
library("ggrepel")
library("dplyr")
library("devtools")
load_all("macrophage-gxe-study/housekeeping/")


#Import GWAS summary list
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt",
                                col_names = c("trait","file_name", "type")) %>%
  dplyr::filter(type == "Autoimmune")

##### import eQTL enrichments #####
#Name all files
file_names = as.list(paste0("results/SL1344/garfield/garfield.perm.",gwas_stats_labeled$file_name, ".out"))
names(file_names) = gwas_stats_labeled$trait

#Import enrichments
enrich_df = purrr::map_df(file_names, ~readr::read_delim(., delim = " "), .id = "trait")
eqtl_filtered_enrichments = dplyr::filter(enrich_df, !(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "RA_2012", "CD_2012", "T1D_2", "T2D_1", "PBC","UC","CD")))

##### import caQTL enrichments #####
#Name all files
file_names = as.list(paste0("results/ATAC/garfield/garfield.perm.",gwas_stats_labeled$file_name, ".out"))
names(file_names) = gwas_stats_labeled$trait

#Import enrichments
enrich_df = purrr::map_df(file_names, ~readr::read_delim(., delim = " "), .id = "trait")
caqtl_filtered_enrichments = dplyr::filter(enrich_df, !(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "RA_2012", "CD_2012", "T1D_2", "T2D_1", "PBC","UC","CD")))

# Combine eQTL and caQTL enrichments
joint_enrichments = dplyr::bind_rows(eqtl_filtered_enrichments, caqtl_filtered_enrichments) %>%
  dplyr::mutate(Pvalue = ifelse(Pvalue == 0, 1e-5, Pvalue)) %>%
  dplyr::mutate(log10p = -log(Pvalue, 10)) %>%
  tidyr::separate(Annotation, into = c("condition_name", "qtl_type"), sep ="_e|_c") %>%
  dplyr::select(-qtl_type) %>%
  dplyr::left_join(figureNames(), by = "condition_name") %>%
  dplyr::mutate(phenotype = ifelse(Type == "eQTL", "RNA-seq", "ATAC-seq"))
joint_enrichments$figure_name = factor(joint_enrichments$figure_name, levels = rev(levels(joint_enrichments$figure_name)))

# Filter by p-value threshold
joint_enrichments_e5 = dplyr::filter(joint_enrichments, PThresh == 1e-5)
joint_enrichments_e8 = dplyr::filter(joint_enrichments, PThresh == 1e-8)

#Calculate max fold enrichments
mean_or_df = dplyr::group_by(joint_enrichments_e5, trait) %>% 
  dplyr::mutate(mean_OR = mean(Beta)) %>% arrange(-mean_OR) %>%
  dplyr::ungroup()
ranked_enrichments = dplyr::mutate(joint_enrichments_e5, trait = factor(trait, levels = unique(mean_or_df$trait))) %>%
  dplyr::mutate(phenotype = factor(phenotype, levels = c("RNA-seq", "ATAC-seq")))

#Make plots
enrichment_plot_e5 = ggplot(ranked_enrichments, aes(x = Beta, xmin = CI95_lower, xmax = CI95_upper, y = figure_name)) + 
  geom_point() + 
  geom_errorbarh(aes(height = 0)) + 
  facet_grid(trait~phenotype) +
  theme_light() +
  ylab("Condition") +
  xlab("log(OR)") +
  geom_vline(aes(xintercept = 0), size = 0.3)
ggsave("figures/supplementary/grafield_enrichment_e-5.pdf", plot = enrichment_plot_e5, width = 4, height = 6)

#Calculate max fold enrichments
mean_or_df = dplyr::group_by(joint_enrichments_e8, trait) %>% 
  dplyr::mutate(mean_OR = mean(Beta)) %>% arrange(-mean_OR) %>%
  dplyr::ungroup()
ranked_enrichments = dplyr::mutate(joint_enrichments_e8, trait = factor(trait, levels = unique(mean_or_df$trait))) %>%
  dplyr::mutate(phenotype = factor(phenotype, levels = c("RNA-seq", "ATAC-seq"))) %>%
  dplyr::filter(trait != "T2D")

#Make plots
enrichment_plot_r8 = ggplot(ranked_enrichments, aes(x = Beta, xmin = CI95_lower, xmax = CI95_upper, y = figure_name)) + 
  geom_point() + 
  geom_errorbarh(aes(height = 0)) + 
  facet_grid(trait~phenotype) +
  theme_light() +
  ylab("Condition") +
  xlab("log(OR)") +
  geom_vline(aes(xintercept = 0), size = 0.3)
ggsave("figures/supplementary/grafield_enrichment_e-8.pdf", plot = enrichment_plot_r8, width = 4, height = 6)

#Export raw enricment p-values and FEs
garfield_table = dplyr::select(joint_enrichments_e5, trait, OR, Pvalue, Beta, CI95_lower, CI95_upper, NAnnotThesh, NAnnot, NThresh, N, condition_name, phenotype)
write.table(garfield_table, "figures/tables/GARFIELD_enrichments.txt", sep = "\t", quote = FALSE, row.names = FALSE)



