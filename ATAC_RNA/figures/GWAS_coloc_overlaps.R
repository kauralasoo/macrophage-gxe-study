library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Functions
identifyColocHits <- function(coloc_df, PP_power_thresh = 0.8, PP_coloc_thresh = 0.9, nsnps_thresh = 10){
  coloc_hits = dplyr::filter(coloc_df, PP_power > PP_power_thresh) %>% 
    dplyr::group_by(trait, gwas_lead, gene_id) %>% 
    dplyr::arrange(trait, gene_id, gwas_lead, -PP.H4.abf) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::filter(PP_coloc > PP_coloc_thresh) %>%
    dplyr::filter(nsnps > nsnps_thresh) %>%
    dplyr::ungroup()
  return(coloc_hits)
}

countConditionSpecificOverlaps <- function(coloc_filtered, PP_power_thresh = 0.8, PP_coloc_thresh = 0.9){
  #Count the number of overlaps added by each additional condition
  coloc_counts = dplyr::mutate(coloc_filtered, is_hit = ifelse(PP_power > PP_power_thresh & PP_coloc > PP_coloc_thresh, 1, 0)) %>% 
    dplyr::select(summarised_trait, gene_id, snp_id, condition_name, is_hit) %>% 
    tidyr::spread(condition_name, is_hit) %>%
    dplyr::mutate(naive_IFNg = pmax(naive, IFNg)) %>%
    dplyr::mutate(naive_IFNg_SL1344 = pmax(naive_IFNg, SL1344)) %>%
    dplyr::mutate(all = pmax(naive_IFNg_SL1344, IFNg_SL1344)) %>%
    dplyr::mutate(IFNg_added = naive_IFNg - naive) %>%
    dplyr::mutate(SL1344_added = naive_IFNg_SL1344 - naive_IFNg) %>%
    dplyr::mutate(IFNg_SL1344_added = all - naive_IFNg_SL1344) %>%
    dplyr::select(summarised_trait, gene_id, snp_id, naive, IFNg_added, SL1344_added, IFNg_SL1344_added) %>%
    dplyr::rename(IFNg = IFNg_added, SL1344 = SL1344_added, IFNg_SL1344 = IFNg_SL1344_added) %>%
    tidyr::gather("condition_name", "is_hit", naive:IFNg_SL1344) %>% 
    dplyr::filter(is_hit == 1) %>%
    dplyr::left_join(figureNames(), by = "condition_name")
}

importAndFilterColocHits <- function(gwas_stats, coloc_suffix = ".eQTL.1e+05.coloc.txt", coloc_prefix = "results/SL1344/coloc/coloc_lists/",
                                     PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, gwas_pval_thresh = 1e-6){
  #Import coloc hits
  
  #Name all files
  file_names = as.list(paste0(coloc_prefix, gwas_stats_labeled$trait, coloc_suffix))
  names(file_names) = gwas_stats_labeled$trait
  
  #Import enrichments
  coloc_df = purrr::map_df(file_names, ~readr::read_delim(., delim = "\t"), .id = "trait") %>%
    dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power) %>%
    dplyr::mutate(summarised_trait = ifelse(trait %in% c("IBD","UC","CD"), "IBD", trait))
  
  #Identify one overlap GWAS lead varaint
  coloc_hits = identifyColocHits(coloc_df, PP_power_thresh, PP_coloc_thresh, nsnps_thresh) %>%
    dplyr::filter(gwas_pval < 1e-6) %>%
    dplyr::filter(!(gene_id %in% mhc_genes$gene_id))
  
  #Merge IBD overlaps together, keep only stronges association per gene
  coloc_hits_merged = dplyr::group_by(coloc_hits, summarised_trait, gene_id) %>% 
    dplyr::arrange(summarised_trait, gene_id, -PP.H4.abf) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup()
  
  #Retreive posterior probabilitites in all conditions
  coloc_filtered = dplyr::semi_join(coloc_df, coloc_hits_merged, by = c("trait", "gene_id", "snp_id"))
  return(list(coloc_filtered = coloc_filtered, coloc_df = coloc_df))
}

#Import expression data
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import atac data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Identify genes in the MHC region that should be excluded
mhc_genes = dplyr::filter(combined_expression_data$gene_metadata, chr == "6", start > 28510120, end < 33480577)
mhc_peaks = dplyr::filter(atac_data$gene_metadata, chr == "6", start > 28510120, end < 33480577)

#Import GWAS traits
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt",
                                     col_names = c("trait","file_name")) %>%
  dplyr::filter(!(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "CD_2012", "RA_2012", "T2D_1", "MS", "T1D", "T1D_2", "PBC")))

#Import unconvincing coloc overlaps that should be filtered out:
unconvincing_coloc = read.table("macrophage-gxe-study/data/gwas_catalog/unconvincing_coloc.txt", stringsAsFactors = FALSE, header = TRUE)

##### eQTL overlaps ####
#Import coloc output
eqtl_100kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".eQTL.1e+05.coloc.txt", 
                                           PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                           gwas_pval_thresh = 1e-6)$coloc_filtered %>%
  dplyr::left_join(gene_name_map, by = "gene_id")
saveRDS(eqtl_100kb_hits, "results/SL1344/coloc/eQTL_coloc_100kb_hits.rds")

eqtl_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".eQTL.2e+05.coloc.txt", 
                                           PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                           gwas_pval_thresh = 1e-6)$coloc_filtered %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>%
  dplyr::anti_join(unconvincing_coloc, by = c("gene_name", "trait"))
saveRDS(eqtl_200kb_hits, "results/SL1344/coloc/eQTL_coloc_200kb_hits.rds")

##### caQTL overlaps ####
#Import coloc output
caqtl_100kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".caQTL.1e+05.coloc.txt", 
                                           PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                           gwas_pval_thresh = 1e-6)$coloc_filtered %>%
  dplyr::mutate(gene_name = gene_id)
saveRDS(caqtl_100kb_hits, "results/SL1344/coloc/caQTL_coloc_100kb_hits.rds")

caqtl_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".caQTL.2e+05.coloc.txt", 
                                           PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                           gwas_pval_thresh = 1e-6)$coloc_filtered %>%
  dplyr::mutate(gene_name = gene_id) %>%
  dplyr::anti_join(unconvincing_coloc, by = c("gene_name", "trait"))
saveRDS(caqtl_200kb_hits, "results/SL1344/coloc/caQTL_coloc_200kb_hits.rds")

#Condense multi-peak caQTLs to a single hit
condensed_caqtl_hits = dplyr::group_by(caqtl_200kb_hits, summarised_trait, gene_id) %>% dplyr::arrange(summarised_trait, -PP.H4.abf) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(summarised_trait, gwas_lead, -PP.H4.abf) %>% 
  dplyr::group_by(summarised_trait, gwas_lead) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()
caqtl_200kb_filtered_hits = dplyr::semi_join(caqtl_200kb_hits, condensed_caqtl_hits, by = c("gene_id", "trait"))
saveRDS(caqtl_200kb_filtered_hits, "results/SL1344/coloc/caQTL_coloc_200kb_hits.rds")


#### Count the number of coloc hits by condition and by trait ####

#Partition into conditions
eqtl_coloc_counts = countConditionSpecificOverlaps(eqtl_200kb_hits, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
eqtl_total_counts = group_by(eqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) %>% 
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "RNA-seq")

#Total counts for caQTLs
caqtl_coloc_counts = countConditionSpecificOverlaps(caqtl_200kb_filtered_hits, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
caqtl_total_counts = dplyr::group_by(caqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) 
caqtl_totals = dplyr::bind_rows(caqtl_total_counts, data_frame(figure_name = c("S","I+S"), overlap_count = c(0,0))) %>%
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "ATAC-seq") %>%
  dplyr::mutate(figure_name = factor(figure_name, levels = c("N","I","S","I+S")))

#Merge both counts
coloc_detection_counts = dplyr::bind_rows(eqtl_total_counts, caqtl_totals)

#Make a barplot with overlap counts
coloc_counts_plot = ggplot(coloc_detection_counts, aes(x = figure_name, y = total_overlap, group = phenotype, color = phenotype)) + 
  geom_point() +
  geom_line() +
  xlab("Condition") + 
  ylab("Cumulative number of overlaps") +
  scale_y_continuous(limits = c(0,25)) +
  theme_light() + 
  scale_color_manual(values = c("#e66101","#5e3c99"), name = "") +
  theme(legend.position = "top")
ggsave("figures/main_figures/coloc_QTL_counts.pdf", plot = coloc_counts_plot, width = 3, height = 3.5)


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
  coord_flip() +
  theme_light() +
  xlab("Trait") +
  ylab("Number of colocalised loci")  + 
  scale_fill_manual(values = c("#e66101","#5e3c99"), name = "", guide = guide_legend(reverse = TRUE)) +
  theme(legend.position = "top") 
ggsave("figures/main_figures/coloc_trait_counts.pdf", plot = coloc_traits_plot, width = 3, height = 3.5)



#Identify and count shared QTLs
shared_qtls = dplyr::semi_join(eqtl_200kb_hits, condensed_caqtl_hits, by = "gwas_lead")

#Count chared eQTL and caQTL overlaps
shared_counts = dplyr::select(shared_qtls, summarised_trait, gene_id) %>% 
  unique() %>% 
  dplyr::group_by(summarised_trait) %>% 
  dplyr::summarise(overlap_count = length(summarised_trait)) %>% 
  dplyr::mutate(phenotype = "shared")


