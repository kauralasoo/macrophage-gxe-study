library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import expression data
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name) %>% dplyr::rename(phenotype_id = gene_id)

#Import atac data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Identify genes in the MHC region that should be excluded
mhc_genes = dplyr::filter(combined_expression_data$gene_metadata, chr == "6", start > 28510120, end < 33480577) %>%
  dplyr::rename(phenotype_id = gene_id)
mhc_peaks = dplyr::filter(atac_data$gene_metadata, chr == "6", start > 28510120, end < 33480577) %>%
  dplyr::rename(phenotype_id = gene_id)

#Import GWAS traits
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt",
                                     col_names = c("trait","file_name")) %>%
  dplyr::filter(!(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "CD_2012", "RA_2012", "T2D_1", "MS", "T1D", "T1D_2", "PBC"))) %>%
  dplyr::filter(trait %in% c("AD","CEL","CD","IBD","UC","NAR","RA","SCZ","SLE","T2D"))

#Import unconvincing coloc overlaps that should be filtered out:
unconvincing_coloc = read.table("macrophage-gxe-study/data/gwas_catalog/unconvincing_coloc.txt", stringsAsFactors = FALSE, header = TRUE)

##### eQTL overlaps ####
#Import coloc output
eqtl_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".SL1344.2e5.txt", 
                                           coloc_prefix = "processed/ATAC_RNA/coloc/",
                                           PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                           gwas_pval_thresh = 1e-6, mhc_phenotypes = mhc_genes)$coloc_filtered %>%
  dplyr::left_join(gene_name_map, by = "phenotype_id") %>%
  dplyr::anti_join(unconvincing_coloc, by = c("gene_name", "trait")) %>%
  dplyr::select(-.row)

#Perform additional coloc 
#Run additional coloc for rescued caQTL-eQTL pairs
rescued_pairs = read.table("results/ATAC_RNA_overlaps/caQTL_eQTL_rescued_pairs.txt",
                           header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::transmute(phenotype_id = gene_id, snp_id, trait)
gwas_stats_prefix = dplyr::mutate(gwas_stats_labeled, 
                                  gwas_prefix = file.path("../../datasets/Inflammatory_GWAS/", file_name))

#Perform coloc
colocs = purrr::map_df(qtlResults()$rna_fastqtl, function(qtl_path){
  res = purrr::by_row(rescued_pairs, ~colocMolecularQTLs(., qtl_summary_path = qtl_path, 
                                                         gwas_summary_path = paste0(dplyr::filter(gwas_stats_prefix, trait == .$trait)$gwas_prefix, ".sorted.txt.gz"), 
                                                         GRCh37_variants = GRCh37_variants, 
                                                         GRCh38_variants = GRCh38_variants, N_qtl = 84, cis_dist = 1e5, QTLTools = FALSE)$summary,.collate = "rows")
  return(res)
}, .id = "condition_name")
rescued_colocs = dplyr::mutate(colocs, PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power) %>%
  dplyr::mutate(summarised_trait = ifelse(trait %in% c("IBD","UC","CD"), "IBD", trait)) 
rescued_hits = identifyColocHits(rescued_colocs) %>% 
  dplyr::anti_join(eqtl_200kb_hits, by = "phenotype_id")
rescued_filtered = dplyr::semi_join(rescued_colocs, rescued_hits, by = c("trait", "phenotype_id", "snp_id")) %>%
  dplyr::left_join(gene_name_map, by = "phenotype_id") %>%
  dplyr::select(-.row)
write.table(rescued_filtered, "figures/tables/eQTL_coloc_100kb_rescued.txt", sep ="\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

saveRDS(eqtl_200kb_hits, "results/SL1344/coloc/eQTL_coloc_200kb_hits.rds")
write.table(eqtl_200kb_hits, "figures/tables/eQTL_coloc_200kb_hits.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)




#Keep one gene per summarised trait
condensed_eQTL_hits = dplyr::group_by(eqtl_200kb_hits, summarised_trait, phenotype_id) %>% 
  dplyr::arrange(summarised_trait, -PP.H4.abf) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(summarised_trait, gwas_lead, -PP.H4.abf) %>% 
  dplyr::group_by(summarised_trait, gwas_lead) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()
eqtl_200kb_filtered_hits = dplyr::semi_join(eqtl_200kb_hits_all, condensed_eQTL_hits, by = c("phenotype_id", "trait")) %>%
  dplyr::filter(gene_name != "FADS2")


##### caQTL overlaps ####
#Import coloc output
caqtl_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".ATAC.2e5.txt", 
                                           coloc_prefix = "processed/ATAC_RNA/coloc/",
                                           PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                           gwas_pval_thresh = 1e-6, mhc_phenotypes = mhc_peaks)$coloc_filtered %>%
  dplyr::mutate(gene_name = phenotype_id) %>%
  dplyr::anti_join(unconvincing_coloc, by = c("gene_name", "trait"))  %>%
  dplyr::select(-.row)

#Condense multi-peak caQTLs to a single hit
condensed_caqtl_hits = dplyr::group_by(caqtl_200kb_hits, summarised_trait, phenotype_id) %>% 
  dplyr::arrange(summarised_trait, -PP.H4.abf) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(summarised_trait, gwas_lead, -PP.H4.abf) %>% 
  dplyr::group_by(summarised_trait, gwas_lead) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()
caqtl_200kb_filtered_hits = dplyr::semi_join(caqtl_200kb_hits, condensed_caqtl_hits, by = c("phenotype_id", "trait"))
saveRDS(caqtl_200kb_filtered_hits, "results/SL1344/coloc/caQTL_coloc_200kb_hits.rds")
write.table(caqtl_200kb_filtered_hits, "figures/tables/caQTL_coloc_200kb_hits.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

#### Count the number of coloc hits by condition and by trait ####

#Partition into conditions
eqtl_coloc_counts = countConditionSpecificOverlaps(eqtl_200kb_filtered_hits, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
eqtl_total_counts = group_by(eqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) %>% 
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "RNA-seq")

#Total counts for caQTLs
caqtl_coloc_counts = countConditionSpecificOverlaps(caqtl_200kb_filtered_hits, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
caqtl_total_counts = dplyr::group_by(caqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) 
caqtl_totals = dplyr::bind_rows(caqtl_total_counts, data_frame(figure_name = c("I+S"), overlap_count = c(0))) %>%
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "ATAC-seq") %>%
  dplyr::mutate(figure_name = factor(figure_name, levels = c("N","I","S","I+S")))

#Merge both counts
coloc_detection_counts = dplyr::bind_rows(eqtl_total_counts, caqtl_totals) %>%
  dplyr::mutate(figure_name = as.character(figure_name)) %>%
  dplyr::mutate(figure_name = ifelse(figure_name == "I", "N\nI", figure_name)) %>%
  dplyr::mutate(figure_name = ifelse(figure_name == "S", "N\nI\nS", figure_name)) %>%
  dplyr::mutate(figure_name = ifelse(figure_name == "I+S", "N\nI\nS\nI+S", figure_name))

#Make a barplot with overlap counts
coloc_counts_plot = ggplot(coloc_detection_counts, aes(x = figure_name, y = total_overlap, group = phenotype, color = phenotype)) + 
  geom_point() +
  geom_line() +
  xlab("Conditions included") + 
  ylab("Cumulative number of overlaps") +
  scale_y_continuous(limits = c(0,25)) +
  theme_light() + 
  scale_color_manual(values = c("#e66101","#5e3c99"), name = "") +
  theme(legend.position = "top")
ggsave("figures/main_figures/coloc_QTL_counts.pdf", plot = coloc_counts_plot, width = 2.6, height = 3)


#Count overlaps by trait
eqtl_by_trait = dplyr::group_by(eqtl_coloc_counts, summarised_trait) %>% 
  dplyr::summarise(overlap_count = length(summarised_trait)) %>% 
  dplyr::mutate(phenotype = "RNA-seq")

caqtl_by_trait = dplyr::group_by(caqtl_coloc_counts, summarised_trait) %>% 
  dplyr::summarise(overlap_count = length(summarised_trait)) %>% 
  dplyr::mutate(phenotype = "ATAC-seq")


#Merge counts
merged_counts = dplyr::bind_rows(eqtl_by_trait, caqtl_by_trait) %>% 
  dplyr::mutate(summarised_trait = factor(summarised_trait, levels = rev(c("IBD","RA","SCZ","CEL","SLE","AD","NAR","T2D"))))

#Make a plot of coloc counts per trait
coloc_traits_plot = ggplot(merged_counts, aes(x = summarised_trait, y = overlap_count, fill = phenotype)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
  coord_flip() +
  theme_light() +
  xlab("Trait") +
  ylab("Number of colocalised loci")  + 
  scale_fill_manual(values = c("#e66101","#5e3c99"), name = "", guide = guide_legend(reverse = TRUE)) +
  theme(legend.position = "top")
ggsave("figures/main_figures/coloc_trait_counts.pdf", plot = coloc_traits_plot, width = 2.6, height = 3)



#Identify and count shared QTLs
shared_qtls = dplyr::semi_join(eqtl_200kb_hits_all, condensed_caqtl_hits, by = "gwas_lead")

#Count chared eQTL and caQTL overlaps
shared_counts = dplyr::select(shared_qtls, summarised_trait, phenotype_id) %>% 
  unique() %>% 
  dplyr::group_by(summarised_trait) %>% 
  dplyr::summarise(overlap_count = length(summarised_trait)) %>% 
  dplyr::mutate(phenotype = "shared")


#Show that ICOSLG eQTL is not colocalised with the UC GWAS hit
a = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".eQTL.2e+05.coloc.txt", 
                         PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                         gwas_pval_thresh = 1e-6)$coloc_df
View(dplyr::filter(a, gene_id == "ENSG00000160223", trait == "UC", snp_id == "rs4819387"))


