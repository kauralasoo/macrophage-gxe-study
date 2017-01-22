library("ggplot2")
library("ggrepel")
library("dplyr")
library("devtools")
load_all("macrophage-gxe-study/housekeeping/")


#Import GWAS summary list
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt",
                                col_names = c("trait","file_name"))

##### import eQTL enrichments #####
#Name all files
file_names = as.list(paste0("results/SL1344/coloc/GARFIELD/eQTL/garfield.perm.",gwas_stats_labeled$file_name, ".out"))
names(file_names) = gwas_stats_labeled$trait

#Import enrichments
enrich_df = purrr::map_df(file_names, ~readr::read_delim(., delim = " "), .id = "trait")
eqtl_filtered_enrichments = dplyr::filter(enrich_df, PThresh == 1e-5) %>%
  dplyr::filter(!(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "RA_2012", "CD_2012", "T1D_2", "T2D_1", "PBC","UC","CD")))

##### import caQTL enrichments #####
#Name all files
file_names = as.list(paste0("results/SL1344/coloc/GARFIELD/caQTL/garfield.perm.",gwas_stats_labeled$file_name, ".out"))
names(file_names) = gwas_stats_labeled$trait

#Import enrichments
enrich_df = purrr::map_df(file_names, ~readr::read_delim(., delim = " "), .id = "trait")
caqtl_filtered_enrichments = dplyr::filter(enrich_df, PThresh == 1e-5) %>%
  dplyr::filter(!(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "RA_2012", "CD_2012", "T1D_2", "T2D_1", "PBC","UC","CD"))) %>%
  dplyr::rename(Annotation = annotation)

# Combine eQTL and caQTL enrichments
joint_enrichments = dplyr::bind_rows(eqtl_filtered_enrichments, caqtl_filtered_enrichments) %>%
  dplyr::mutate(EmpPval = ifelse(EmpPval == 0, 1e-5, EmpPval)) %>%
  dplyr::mutate(log10p = -log(EmpPval, 10)) %>%
  tidyr::separate(Annotation, into = c("condition_name", "qtl_type"), sep ="_e|_c") %>%
  dplyr::select(-qtl_type) %>%
  dplyr::left_join(figureNames(), by = "condition_name") %>%
  dplyr::mutate(phenotype = ifelse(Type == "eQTL", "RNA-seq", "ATAC-seq"))

#Calculate max fold enrichments
mean_fe_df = dplyr::group_by(joint_enrichments, trait) %>% 
  dplyr::mutate(mean_FE = mean(FE)) %>% arrange(mean_FE) %>%
  dplyr::ungroup()
ranked_enrichments = dplyr::mutate(joint_enrichments, trait = factor(trait, levels = unique(mean_fe_df$trait))) %>%
  dplyr::mutate(phenotype = factor(phenotype, levels = c("RNA-seq", "ATAC-seq")))

#Make plots
enrichment_plot = ggplot(ranked_enrichments, aes(x = FE, y = trait, color = figure_name, size = log10p)) + 
  geom_jitter(height = 0, width = 2) +
  scale_x_continuous(limits = c(0, 41)) + 
  theme_light() +
  scale_size(range = c(1,4), name = "-log10\np-value") +
  facet_wrap(~phenotype) +
  xlab("Fold enrichment") +
  ylab("Trait") +
  scale_color_manual(values = conditionPalette(), name="Condition") + 
  theme(strip.text.x = element_text(colour = "grey10"), strip.background = element_rect(fill = "grey85"))
  
ggsave("figures/main_figures/grafield_enrichment.pdf", plot = enrichment_plot, width = 4.5, height = 3.5)

#Export raw enricment p-values and FEs
garfield_table = dplyr::select(joint_enrichments, trait, FE, EmpPval, NAnnotThesh, NAnnot, NThresh, N, condition_name, phenotype)
write.table(garfield_table, "figures/tables/GARFIELD_enrichments.txt", sep = "\t", quote = FALSE, row.names = FALSE)



