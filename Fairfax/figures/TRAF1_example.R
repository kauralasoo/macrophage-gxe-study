library("dplyr")
library("tidyr")
library("purrr")
library("UpSetR")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("../wiggleplotr/")

#Import SummarizedExperiment
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")
expression_mat = assays(se_fairfax)$exprs
sample_metadata = colData(se_fairfax) %>% tbl_df2()
gene_metadata = rowData(se_fairfax) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id)

#Import genotypes
vcf_file = seqUtils::gdsToMatrix("processed/Fairfax/geno_by_chr/9.gds")
variant_information = importVariantInformation("processed/Fairfax/merged_genotypes/fairfax_genotypes.variant_information.txt.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name"))

#Make a QTL plot
TRAF1_data = constructQtlPlotDataFrame("ILMN_1698218", "rs10985070", 
                                       expression_mat, 
                                       vcf_file$genotypes, 
                                       sample_metadata, 
                                       gene_metadata) %>%
  dplyr::left_join(constructGenotypeText("rs10985070", variant_information), by = "genotype_value")

TRAF1_qtl_plot = plotQtlCol(TRAF1_data, scales = "free_y")
ggsave("figures/supplementary/fairfax_TRAF1_boxplot.pdf", width = 4, height = 6)



#Make a manhattan plot of the region

#Import QTL and GWAS summary stats and convert them to the same GRCh38 coordinate space
qtl_df = data_frame(phenotype_id = "ILMN_1698218", snp_id = "rs10985070", trait = "RA")
qtl_paths = list(CD14 = "/Volumes/JetDrive/databases/Fairfax/full/CD14.nominal.sorted.txt.gz",
                 IFN = "/Volumes/JetDrive/databases/Fairfax/full/IFN.nominal.sorted.txt.gz",
                 LPS2 = "/Volumes/JetDrive/databases/Fairfax/full/LPS2.nominal.sorted.txt.gz",
                 LPS24 = "/Volumes/JetDrive/databases/Fairfax/full/LPS24.nominal.sorted.txt.gz")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, gwas_dir = "/Volumes/JetDrive/datasets/Inflammatory_GWAS/", qtl_paths = qtl_paths, 
                                         GRCh37_variants = variant_information, GRCh38_variants = variant_information, cis_dist = 2e5, type = "QTLTools") %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::rename(track_id = condition_name) %>%
  dplyr::mutate(track_id = factor(track_id, levels = c("RA","CD14","IFN","LPS2","LPS24")))

#Make a manhattan plot
region_coords = c(min(qtl_summary$pos), max(qtl_summary$pos))
gwas_plot = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id == "RA"), limits = region_coords, color_R2 = TRUE)
qtl_plot = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id != "RA"), limits = region_coords, color_R2 = TRUE, data_track = FALSE) +
  theme(plot.margin=unit(c(0.1,1,0.1,1),"line"),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(colour = "grey10"),
        strip.background = element_rect(fill = "grey85")) +
  xlab("Chromosome 9 position")

joint_plot = cowplot::plot_grid(gwas_plot, qtl_plot, 
                                align = "v", ncol = 1, rel_heights = c(1,4.5))
ggsave("figures/supplementary//fairfax_TRAF1_manhattan.pdf", joint_plot, width = 4, height = 6)


