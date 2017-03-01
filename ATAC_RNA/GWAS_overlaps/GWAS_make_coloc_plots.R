library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("ggplot2")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Functions
plotColoc <- function(df, plot_title){
  plot = ggplot(df, aes(x = pos, y = log10p)) + 
    geom_point() + 
    facet_grid(condition_name~.) +
    labs(title = plot_title) + 
    theme_light()
  return(plot)
}

#Import expression data
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import old and new variant coordinates
GRCh38_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name"))

#### Import list of colocalised eQTLs in 200kb window ####
eqtl_coloc_hits_200kb = readRDS("results/SL1344/coloc/eQTL_coloc_200kb_hits.rds") %>%
  dplyr::select(phenotype_id, snp_id, trait, gwas_lead, gene_name) %>% unique() %>%
  dplyr::mutate(plot_title = paste(trait, gene_name, gwas_lead, sep = "_"))

#Fetch data for all eqtl hits
plot_data = purrr::by_row(eqtl_coloc_hits_200kb, 
                          ~importSummariesForPlotting(., gwas_stats_labeled, gwas_dir = "~/datasets/Inflammatory_GWAS/",
                                                      qtl_paths = qtlResults()$rna_fastqtl, 
                                                      GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants, cis_dist = 2e5), .to = "data")

#Make plots
plots = purrr::by_row(plot_data, ~plotColoc(.$data[[1]], .$plot_title), .to = "plot")
plot_list = plots$plot
names(plot_list) = plots$plot_title
savePlotList(plot_list, "processed/ATAC_RNA/coloc_plots/eQTLs/")


#Import list of colocalised caQTLs in 200kb window
caqtl_coloc_hits_200kb = readRDS("results/SL1344/coloc/caQTL_coloc_200kb_hits.rds") %>%
  dplyr::select(phenotype_id, snp_id, trait, gwas_lead, gene_name) %>% unique() %>%
  dplyr::mutate(plot_title = paste(trait, gene_name, gwas_lead, sep = "_"))

#Fetch data for all eqtl hits
plot_data = purrr::by_row(caqtl_coloc_hits_200kb, 
                          ~importSummariesForPlotting(., gwas_stats_labeled, gwas_dir = "~/datasets/Inflammatory_GWAS/",
                                                      qtl_paths = qtlResults()$atac_fastqtl, 
                        GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants, cis_dist = 2e5), .to = "data")

#Make plots
plots = purrr::by_row(plot_data, ~plotColoc(.$data[[1]], .$plot_title), .to = "plot")
plot_list = plots$plot
names(plot_list) = plots$plot_title
savePlotList(plot_list, "processed/ATAC_RNA/coloc_plots/caQTLs/")

