library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import SummarizedExperiments
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
se_reviseAnnotations = readRDS("results/acLDL/acLDL_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/acLDL/acLDL_leafcutter_counts.rds")
se_featureCounts = readRDS("results/acLDL/acLDL_combined_expression_data_covariates.se.rds")

#Identify genes in the MHC region that should be excluded
mhc_ensembl = dplyr::filter(tbl_df2(rowData(se_ensembl)), chr == "6", transcript_start > 28510120, transcript_start < 33480577) %>%
  dplyr::rename(phenotype_id = transcript_id)
mhc_revised = dplyr::filter(tbl_df2(rowData(se_reviseAnnotations)), chr == "6", start > 28510120, end < 33480577) %>%
  dplyr::rename(phenotype_id = transcript_id)
mhc_leafcutter = dplyr::filter(tbl_df2(rowData(se_leafcutter)), chr == "6", start > 28510120, end < 33480577) %>%
  dplyr::rename(phenotype_id = transcript_id)
mhc_featureCounts = dplyr::filter(tbl_df2(rowData(se_featureCounts)), chr == "6", start > 28510120, end < 33480577) %>%
  dplyr::rename(phenotype_id = gene_id)


#Gene names
ensembl_name_map = dplyr::select(tbl_df2(rowData(se_ensembl)), transcript_id, gene_name) %>% dplyr::rename(phenotype_id = transcript_id)
revised_name_map = dplyr::select(tbl_df2(rowData(se_reviseAnnotations)), transcript_id, gene_name) %>% dplyr::rename(phenotype_id = transcript_id)
leafcutter_name_map = dplyr::select(tbl_df2(rowData(se_leafcutter)), transcript_id, gene_name) %>% dplyr::rename(phenotype_id = transcript_id)
featureCounts_name_map = dplyr::select(tbl_df2(rowData(se_featureCounts)), gene_id, gene_name) %>% dplyr::rename(phenotype_id = gene_id)

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

#featureCounts
featureCounts_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".featureCounts.2e5.txt", 
                                                 coloc_prefix = "processed/acLDL/coloc/",
                                                 PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, 
                                                 gwas_pval_thresh = 1e-6, mhc_phenotypes = mhc_featureCounts)$coloc_filtered %>%
  dplyr::left_join(featureCounts_name_map, by = "phenotype_id") %>%
  #dplyr::anti_join(unconvincing_coloc, by = c("gene_name", "trait")) %>%
  dplyr::select(-.row)

#Put all GWAS overlaps into a single list
gwas_olaps = list(ensembl_87 = ensembl_200kb_hits, revisedAnnotation = revised_200kb_hits, 
                  leafcutter = leafcutter_200kb_hits, featureCounts = featureCounts_200kb_hits)
saveRDS(gwas_olaps, "acLDL_figures/tables/GWAS_coloc_hits.rds")

#Save a text table as well
gwas_olaps_tbl = purrr::map_df(gwas_olaps, identity ,.id = "annotation")
write.table(gwas_olaps_tbl, "acLDL_figures/tables/GWAS_coloc_hits.txt", sep = "\t", quote = FALSE, row.names = FALSE)






###### Make plots for all of the genes ######
#Import variant information
GRCh38_variants = importVariantInformation("genotypes/acLDL/imputed_20151005/imputed.70_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/acLDL/imputed_20151005/GRCh37/imputed.70_samples.variant_information.GRCh37.txt.gz")

#Filter coloc hits
ensembl_hits = dplyr::select(ensembl_200kb_hits, phenotype_id, snp_id, trait, gwas_lead, gene_name) %>% 
  unique() %>%
  dplyr::mutate(plot_title = paste(trait, gene_name, gwas_lead, sep = "_"))

#Define location of summary files
ensembl_summaries = list(Ctrl = "processed/acLDL/fastqtl_output/ensembl_87/sorted/Ctrl.nominal.sorted.txt.gz",
     AcLDL = "processed/acLDL/fastqtl_output/ensembl_87/sorted/AcLDL.nominal.sorted.txt.gz")

#Fetch data for all eqtl hits
plot_data = purrr::by_row(ensembl_hits, 
                          ~importSummariesForPlotting(., gwas_stats_labeled, 
                                                      gwas_dir = "~/datasets/Inflammatory_GWAS/",
                                                      qtl_paths = ensembl_summaries, 
                                                      GRCh37_variants = GRCh37_variants, 
                                                      GRCh38_variants = GRCh38_variants, 
                                                      cis_dist = 2e5, type = "QTLTools"), .to = "data")

#Make plots
plots = purrr::by_row(plot_data, ~plotColoc(.$data[[1]], .$plot_title), .to = "plot")
plot_list = setNames(plots$plot, plots$plot_title)
savePlotList(plot_list, "processed/acLDL/coloc_plots/ensembl_87/")


#Filter reviseAnnotations coloc hits
revised_hits = dplyr::select(revised_200kb_hits, phenotype_id, snp_id, trait, gwas_lead, gene_name) %>% 
  unique() %>%
  dplyr::mutate(plot_title = paste(trait, gene_name, gwas_lead, sep = "_"))

#Define location of summary files
revised_summaries = list(Ctrl = "processed/acLDL/fastqtl_output/reviseAnnotations/sorted/Ctrl.nominal.sorted.txt.gz",
                         AcLDL = "processed/acLDL/fastqtl_output/reviseAnnotations/sorted/AcLDL.nominal.sorted.txt.gz")

#Fetch data for all eqtl hits
plot_data = purrr::by_row(revised_hits, 
                          ~importSummariesForPlotting(., gwas_stats_labeled, 
                                                      gwas_dir = "~/datasets/Inflammatory_GWAS/",
                                                      qtl_paths = revised_summaries, 
                                                      GRCh37_variants = GRCh37_variants, 
                                                      GRCh38_variants = GRCh38_variants, 
                                                      cis_dist = 2e5, type = "QTLTools"), .to = "data")

#Make plots
plots = purrr::by_row(plot_data, ~plotColoc(.$data[[1]], .$plot_title), .to = "plot")
plot_list = setNames(plots$plot, plots$plot_title)
savePlotList(plot_list, "processed/acLDL/coloc_plots/reviseAnnotations/")

#Filter leafcutter hits
leafcutter_hits = dplyr::transmute(leafcutter_200kb_hits, phenotype_id, snp_id, trait, gwas_lead, gene_name) %>% 
  unique() %>%
  dplyr::mutate(plot_title = paste(trait, gene_name, gwas_lead, sep = "_"))

#Define location of summary files
leadcutter_summaries = list(Ctrl = "processed/acLDL/fastqtl_output/leafcutter/sorted/Ctrl.nominal.sorted.txt.gz",
                         AcLDL = "processed/acLDL/fastqtl_output/leafcutter/sorted/AcLDL.nominal.sorted.txt.gz")

#Fetch data for all eqtl hits
plot_data = purrr::by_row(leafcutter_hits, 
                          ~importSummariesForPlotting(., gwas_stats_labeled, 
                                                      gwas_dir = "~/datasets/Inflammatory_GWAS/",
                                                      qtl_paths = leadcutter_summaries, 
                                                      GRCh37_variants = GRCh37_variants, 
                                                      GRCh38_variants = GRCh38_variants, 
                                                      cis_dist = 2e5, type = "QTLTools"), .to = "data")

#Make plots
plots = purrr::by_row(plot_data, ~plotColoc(.$data[[1]], .$plot_title), .to = "plot")
plots = dplyr::mutate(plots, plot_title = paste(trait, gene_name, gwas_lead, sep = "_"))
plot_list = setNames(plots$plot, plots$plot_title)
savePlotList(plot_list, "processed/acLDL/coloc_plots/leafcutter/")


#Filter coloc hits
featureCounts_hits = dplyr::select(featureCounts_200kb_hits, phenotype_id, snp_id, trait, gwas_lead, gene_name) %>% 
  unique() %>%
  dplyr::mutate(plot_title = paste(trait, gene_name, gwas_lead, sep = "_"))

#Define location of summary files
featureCounts_summaries = list(Ctrl = "processed/acLDL/fastqtl_output/featureCounts//sorted/Ctrl.nominal.sorted.txt.gz",
                         AcLDL = "processed/acLDL/fastqtl_output/featureCounts/sorted/AcLDL.nominal.sorted.txt.gz")

#Fetch data for all eqtl hits
plot_data = purrr::by_row(featureCounts_hits, 
                          ~importSummariesForPlotting(., gwas_stats_labeled, 
                                                      gwas_dir = "~/datasets/Inflammatory_GWAS/",
                                                      qtl_paths = featureCounts_summaries, 
                                                      GRCh37_variants = GRCh37_variants, 
                                                      GRCh38_variants = GRCh38_variants, 
                                                      cis_dist = 2e5, type = "QTLTools"), .to = "data")

#Make plots
plots = purrr::by_row(plot_data, ~plotColoc(.$data[[1]], .$plot_title), .to = "plot")
plot_list = setNames(plots$plot, plots$plot_title)
savePlotList(plot_list, "processed/acLDL/coloc_plots/featureCounts/")




