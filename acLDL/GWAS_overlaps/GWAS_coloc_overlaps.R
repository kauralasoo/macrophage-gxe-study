library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import SummarizedExperiments
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
se_reviseAnnotations = readRDS("results/acLDL/acLDL_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/acLDL/acLDL_leafcutter_counts.rds")

#Identify genes in the MHC region that should be excluded
mhc_ensembl = dplyr::filter(tbl_df2(rowData(se_ensembl)), chr == "6", transcript_start > 28510120, transcript_start < 33480577) %>%
  dplyr::rename(phenotype_id = transcript_id)
mhc_revised = dplyr::filter(tbl_df2(rowData(se_reviseAnnotations)), chr == "6", start > 28510120, end < 33480577) %>%
  dplyr::rename(phenotype_id = transcript_id)
mhc_leafcutter = dplyr::filter(tbl_df2(rowData(se_leafcutter)), chr == "6", start > 28510120, end < 33480577) %>%
  dplyr::rename(phenotype_id = transcript_id)

#Gene names
ensembl_name_map = dplyr::select(tbl_df2(rowData(se_ensembl)), transcript_id, gene_name) %>% dplyr::rename(phenotype_id = transcript_id)
revised_name_map = dplyr::select(tbl_df2(rowData(se_reviseAnnotations)), transcript_id, gene_name) %>% dplyr::rename(phenotype_id = transcript_id)
leafcutter_name_map = dplyr::select(tbl_df2(rowData(se_leafcutter)), transcript_id, ensembl_gene_id) %>% 
  dplyr::left_join(dplyr::transmute(tbl_df2(rowData(se_ensembl)), ensembl_gene_id = gene_id, gene_name), by = "ensembl_gene_id") %>%
  dplyr::rename(phenotype_id = transcript_id) %>%
  unique()


#Import GWAS traits
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt",
                                     col_names = c("trait","file_name")) %>%
  dplyr::filter(!(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "CD_2012", "RA_2012", "T2D_1", "MS", "T1D", "T1D_2", "PBC")))

#Import coloc output
#Ensembl
ensembl_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".ensembl_87.2e5.txt", 
                                           coloc_prefix = "processed/acLDL/coloc/",
                                           PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                           gwas_pval_thresh = 1e-6, mhc_phenotypes = mhc_ensembl)$coloc_filtered %>%
  dplyr::left_join(ensembl_name_map, by = "phenotype_id") %>%
  #dplyr::anti_join(unconvincing_coloc, by = c("gene_name", "trait")) %>%
  dplyr::select(-.row)

#revisedAnnotations
revised_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".reviseAnnotations.2e5.txt", 
                                              coloc_prefix = "processed/acLDL/coloc/",
                                              PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                              gwas_pval_thresh = 1e-6, mhc_phenotypes = mhc_revised)$coloc_filtered %>%
  dplyr::left_join(revised_name_map, by = "phenotype_id") %>%
  #dplyr::anti_join(unconvincing_coloc, by = c("gene_name", "trait")) %>%
  dplyr::select(-.row)

#revisedAnnotations
leafcutter_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".leafcutter.2e5.txt", 
                                              coloc_prefix = "processed/acLDL/coloc/",
                                              PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                              gwas_pval_thresh = 1e-6, mhc_phenotypes = mhc_leafcutter)$coloc_filtered %>%
  dplyr::left_join(leafcutter_name_map, by = "phenotype_id") %>%
  #dplyr::anti_join(unconvincing_coloc, by = c("gene_name", "trait")) %>%
  dplyr::select(-.row)



