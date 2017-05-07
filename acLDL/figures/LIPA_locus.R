library("dplyr")
library("ggplot2")
library("devtools")
library("rtracklayer")
library("wiggleplotr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../reviseAnnotations/")
load_all("../seqUtils/")

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


#Extract QTL
leafcutter_qtl = dplyr::filter(trqtl_min_pvalues$leafcutter$Ctrl, group_id == "clu_20804")
leafcutter_tx_names = dplyr::filter(gene_meta_list$leafcutter, gene_id == "clu_20804")$transcript_id
leafcutter_tx = leafcutter_granges[leafcutter_tx_names]

#Extract annotatated transcripts
lipa_tx = exons["ENST00000336233"] %>% removeMetadata()
lipa_cdss = cdss["ENST00000336233"] %>% removeMetadata()

#Make a coverage plot
track_data = wiggleplotrGenotypeColourGroup(str1_df, leafcutter_qtl$snp_id, vcf_file$genotypes, 1)
lipa_coverage = wiggleplotr::plotCoverage(c(lipa_tx, leafcutter_tx), lipa_cdss, track_data = track_data, 
                          fill_palette = getGenotypePalette(), 
                          plot_fraction = 0.2, 
                          coverage_type = "line", 
                          rescale_introns = TRUE)
ggsave("acLDL_figures/LIPA_leafcutter_coverage.pdf", lipa_coverage, width = 10, height = 6)

#Make coverage plot in the original scale
lipa_coverage = wiggleplotr::plotCoverage(c(lipa_tx, leafcutter_tx), lipa_cdss, track_data = track_data, 
                                          fill_palette = getGenotypePalette(), 
                                          plot_fraction = 0.2, 
                                          coverage_type = "line", 
                                          rescale_introns = FALSE)
ggsave("acLDL_figures/LIPA_leafcutter_coverage.unscaled.pdf", lipa_coverage, width = 10, height = 6)

#Make a QTL boxplot
sample_metadata = sample_meta_list$leafcutter %>% dplyr::mutate(condition_name = factor(condition_name, levels = c("Ctrl", "AcLDL")))
data = constructQtlPlotDataFrame("10:89247649:89251707:clu_20804", "rs1332328", assays(se_leafcutter)$tpm_ratios, vcf_file$genotypes, 
                                 sample_metadata, 
                                 tbl_df2(rowData(se_leafcutter)) %>% dplyr::mutate(gene_id = transcript_id, gene_name = transcript_id)) %>%
  dplyr::left_join(constructGenotypeText("rs1332328", variant_information), by = "genotype_value")
boxplot = plotQtlRow(data) + 
  ylab("Relative expression")
ggsave("acLDL_figures/LIPA_sQTL_boxplot.pdf", boxplot, width = 4, height = 3)

