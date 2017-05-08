library("dplyr")
library("ggplot2")
library("devtools")
library("rtracklayer")
library("wiggleplotr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../reviseAnnotations/")
load_all("../seqUtils/")
load_all("../wiggleplotr/")

#Import gene exrression data
acldl_list = readRDS("results/acLDL/acLDL_combined_expression_data_covariates.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.rds")
variant_information = importVariantInformation("genotypes/acLDL/imputed_20151005/imputed.70_samples.variant_information.txt.gz")

#Import QTL mapping results
trqtl_min_pvalues = readRDS("results/acLDL/trQTLs/trQTL_min_pvalues.other_tx.rds")

#Import SummarizedExperiments for all phenotypes
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
se_revised = readRDS("results/acLDL/acLDL_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/acLDL/acLDL_leafcutter_counts.rds")
se_list = list(ensembl_87 = se_ensembl, revisedAnnotation = se_revised, leafcutter = se_leafcutter)

#Extract sample metadata
sample_meta_list = purrr::map(se_list, ~colData(.) %>% tbl_df2())
gene_meta_list = purrr::map(se_list, ~rowData(.) %>% tbl_df2())

#Import Ensembl transcript annotations
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

#Leafcutter txs
leafcutter_granges = readRDS("results/acLDL/acLDL_leafcutter.GRangesList.rds")
revised_granges = readRDS("results/reviseAnnotations/reviseAnnotations.GRangesList.rds") %>%
  purrr::flatten()

#Set up sample coverage df
str1_df = wiggleplotrConstructMetadata(acldl_list$counts, 
                                       sample_meta_list$revisedAnnotation, 
                                       "/Volumes/JetDrive/bigwig/", 
                                       bigWig_suffix = ".str1.bw",
                                       condition_name_levels = c("Ctrl","AcLDL"))
str2_df = wiggleplotrConstructMetadata(acldl_list$counts, 
                                       sample_meta_list$revisedAnnotation, 
                                       "/Volumes/JetDrive/bigwig/", 
                                       bigWig_suffix = ".str2.bw",
                                       condition_name_levels = c("Ctrl","AcLDL"))

#### Visualise GWAS overlaps ####
gwas_olaps = readRDS("acLDL_figures/tables/GWAS_coloc_hits.rds")

#CD33 example
cd33_qtl = dplyr::filter(trqtl_min_pvalues$revisedAnnotations$Ctrl, group_id == "ENSG00000105383.contained")
cd33_tx_names = dplyr::filter(gene_meta_list$revisedAnnotation, gene_id == "ENSG00000105383.contained")$transcript_id
cd33_tx = revised_granges[cd33_tx_names]

#Remove the last two exons from both tx
cd33_tx[[1]] = (cd33_tx[[1]])[1:4]
cd33_tx[[2]] = (cd33_tx[[2]])[1:3]
region_coords = c(51224000,51227000)

#Make a coverage plot
track_data = wiggleplotrGenotypeColourGroup(str2_df, cd33_qtl$snp_id, vcf_file$genotypes, 1)
cd33_coverage = wiggleplotr::plotCoverage(cd33_tx, track_data = track_data, 
                                          fill_palette = getGenotypePalette(), 
                                          plot_fraction = 0.2, 
                                          coverage_type = "line", 
                                          rescale_introns = FALSE, region_coords = region_coords, return_subplots_list = TRUE)

#QTL summary list
qtl_summary_list = list(Ctrl = "processed/acLDL/fastqtl_output/reviseAnnotations/sorted/Ctrl.nominal.sorted.txt.gz",
                        AcLDL = "processed/acLDL/fastqtl_output/reviseAnnotations/sorted/AcLDL.nominal.sorted.txt.gz")

#Import gene p-values from the region
qtl_ranges = constructVariantRanges(cd33_qtl, variant_information, cis_dist = 500000)
gene_pvalues = purrr::map_df(qtl_summary_list, ~qtltoolsTabixFetchPhenotypes(qtl_ranges,
                                                                             tabix_file = .)[[1]],
                             .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("Ctrl","AcLDL"))) %>%
  dplyr::mutate(track_id = condition_name) %>%
  dplyr::mutate(pos = snp_start) %>%
  dplyr::arrange(p_nominal) %>%
  addR2FromLead(vcf_file$genotypes)

#Make Manhattan plots
gene_manhattan = makeManhattanPlot(gene_pvalues, limits = region_coords, color_R2 = TRUE, data_track = TRUE)

joint_plot = cowplot::plot_grid(gene_manhattan, cd33_coverage$coverage_plot, cd33_coverage$tx_structure, 
                                align = "v", ncol = 1, rel_heights = c(3,3,2))
ggsave("acLDL_figures/CD33_coverage.pdf", plot = joint_plot, width = 5, height = 6)

#Make a QTL boxplot
sample_metadata = sample_meta_list$revisedAnnotation %>% dplyr::mutate(condition_name = factor(condition_name, levels = c("Ctrl", "AcLDL")))
data = constructQtlPlotDataFrame("ENSG00000105383.clique_1.contained.ENST00000436584", "rs12459419", assays(se_revised)$tpm_ratios, vcf_file$genotypes, 
                                 sample_metadata, 
                                 tbl_df2(rowData(se_revised)) %>% dplyr::mutate(gene_id = transcript_id, gene_name = transcript_id)) %>%
  dplyr::left_join(constructGenotypeText("rs12459419", variant_information), by = "genotype_value")
boxplot = plotQtlRow(data)
ggsave("acLDL_figures/CD33_boxplot.pdf", plot = boxplot, width = 5, height = 6)







#BIN1 example
cd33_qtl = dplyr::filter(trqtl_min_pvalues$revisedAnnotations$Ctrl, group_id == "ENSG00000136717.contained")
cd33_tx_names = dplyr::filter(gene_meta_list$ensembl_87, gene_id == "ENSG00000136717")$transcript_id
cd33_tx = exons[cd33_tx_names]

region_coords = c(127047000,127150000)

#Make a coverage plot
track_data = wiggleplotrGenotypeColourGroup(str2_df, cd33_qtl$snp_id, vcf_file$genotypes, 1)
cd33_coverage = wiggleplotr::plotCoverage(cd33_tx, track_data = track_data, 
                                          fill_palette = getGenotypePalette(), 
                                          plot_fraction = 0.2, 
                                          coverage_type = "line", 
                                          rescale_introns = FALSE, region_coords = region_coords, return_subplots_list = TRUE)

#QTL summary list
qtl_summary_list = list(Ctrl = "processed/acLDL/fastqtl_output/reviseAnnotations/sorted/Ctrl.nominal.sorted.txt.gz",
                        AcLDL = "processed/acLDL/fastqtl_output/reviseAnnotations/sorted/AcLDL.nominal.sorted.txt.gz")

#Import gene p-values from the region
qtl_ranges = constructVariantRanges(cd33_qtl, variant_information, cis_dist = 500000)
gene_pvalues = purrr::map_df(qtl_summary_list, ~qtltoolsTabixFetchPhenotypes(qtl_ranges,
                                                                             tabix_file = .)[[1]],
                             .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("Ctrl","AcLDL"))) %>%
  dplyr::mutate(track_id = condition_name) %>%
  dplyr::mutate(pos = snp_start) %>%
  dplyr::arrange(p_nominal) %>%
  addR2FromLead(vcf_file$genotypes)

#Make Manhattan plots
gene_manhattan = makeManhattanPlot(gene_pvalues, limits = region_coords, color_R2 = TRUE, data_track = TRUE)

joint_plot = cowplot::plot_grid(gene_manhattan, cd33_coverage$coverage_plot, cd33_coverage$tx_structure, 
                                align = "v", ncol = 1, rel_heights = c(3,3,4))
ggsave("acLDL_figures/BIN1_coverage.pdf", plot = joint_plot, width = 10, height = 10)


