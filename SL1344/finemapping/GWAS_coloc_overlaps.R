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
    dplyr::select(trait, gene_id, snp_id, condition_name, is_hit) %>% 
    tidyr::spread(condition_name, is_hit) %>%
    dplyr::mutate(naive_IFNg = pmax(naive, IFNg)) %>%
    dplyr::mutate(naive_IFNg_SL1344 = pmax(naive_IFNg, SL1344)) %>%
    dplyr::mutate(all = pmax(naive_IFNg_SL1344, IFNg_SL1344)) %>%
    dplyr::mutate(IFNg_added = naive_IFNg - naive) %>%
    dplyr::mutate(SL1344_added = naive_IFNg_SL1344 - naive_IFNg) %>%
    dplyr::mutate(IFNg_SL1344_added = all - naive_IFNg_SL1344) %>%
    dplyr::select(trait, gene_id, snp_id, naive, IFNg_added, SL1344_added, IFNg_SL1344_added) %>%
    dplyr::rename(IFNg = IFNg_added, SL1344 = SL1344_added, IFNg_SL1344 = IFNg_SL1344_added) %>%
    tidyr::gather("condition_name", "is_hit", naive:IFNg_SL1344) %>% dplyr::filter(is_hit == 1) %>%
    dplyr::left_join(figureNames(), by = "condition_name")
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
  dplyr::filter(!(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "RA_2012", "CD_2012", "T2D_1", "MS", "T1D", "T1D_2")))# %>%
  #dplyr::filter(!(trait %in% c("UC","CD", "NAR"))) #Remove UC and CD, because there might be a lot of sharing with IBD

##### eQTL overlaps ####
#Import coloc output
#Name all files
file_names = as.list(paste0("results/SL1344/coloc/coloc_lists/", gwas_stats_labeled$trait, ".coloc.txt"))
names(file_names) = gwas_stats_labeled$trait

#Import enrichments
coloc_df = purrr::map_df(file_names, ~readr::read_delim(., delim = "\t"), .id = "trait") %>%
  dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power)

#Identify one overlap GWAS lead varaint
coloc_hits = identifyColocHits(coloc_df, PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 10) %>%
  dplyr::filter(gwas_pval < 1e-6) %>%
  dplyr::filter(!(gene_id %in% mhc_genes$gene_id))

#Retreive posterior probabilitites in all conditions
coloc_filtered = dplyr::semi_join(coloc_df, coloc_hits, by = c("trait", "gene_id", "snp_id")) %>%
  dplyr::left_join(gene_name_map, by = "gene_id")
saveRDS(coloc_filtered, "results/SL1344/coloc/eQTL_coloc_posterior_hits.rds")

#Partition into conditions
coloc_counts = countConditionSpecificOverlaps(coloc_filtered, PP_power_thresh = 0.8, PP_coloc_thresh = .9)

#Make a barplot with overlap counts
eqtl_coloc_counts = ggplot(coloc_counts, aes(x = figure_name, fill = trait)) + 
  geom_bar() +
  xlab("Condition")
ggsave("figures/main_figures/coloc_eQTL_counts.pdf", plot = eqtl_coloc_counts, width = 4.5, height = 4)



##### caQTL overlaps #####
#Import coloc output
#Name all files
file_names = as.list(paste0("results/SL1344/coloc/coloc_lists/", gwas_stats_labeled$trait, ".caQTL.coloc.txt"))
names(file_names) = gwas_stats_labeled$trait

#Import enrichments
caqtl_coloc_df = purrr::map_df(file_names, ~readr::read_delim(., delim = "\t"), .id = "trait") %>%
  dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power)

#Identify one overlap GWAS lead varaint
caqtl_coloc_hits = identifyColocHits(caqtl_coloc_df, PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 10) %>%
  dplyr::filter(gwas_pval < 1e-6) %>%
  dplyr::filter(!(gene_id %in% mhc_peaks$gene_id))

#Retreive posterior probabilitites in all conditions
caqtl_coloc_filtered = dplyr::semi_join(caqtl_coloc_df, caqtl_coloc_hits, by = c("trait", "gene_id", "snp_id")) %>%
  dplyr::mutate(gene_name = gene_id)
saveRDS(caqtl_coloc_filtered, "results/SL1344/coloc/caQTL_coloc_posterior_hits.rds")

#Partition into conditions
caqtl_coloc_counts = countConditionSpecificOverlaps(caqtl_coloc_filtered, PP_power_thresh = 0.8, PP_coloc_thresh = .9)

#Make a barplot with overlap counts
caqtl_coloc_plot = ggplot(caqtl_coloc_counts, aes(x = figure_name, fill = trait)) + 
  geom_bar() +
  xlab("Condition")
ggsave("figures/main_figures/coloc_caQTL_counts.pdf", plot = caqtl_coloc_plot, width = 4.5, height = 4)
