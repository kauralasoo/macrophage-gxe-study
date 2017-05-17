library("dplyr")
library("data.table")
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
gwas_olaps = readRDS("acLDL_figures/tables/GWAS_coloc_hits.rds") %>% 
  purrr::map_df(identity, .id = "phenotype")

#Find FCGR2A overlaps
hits = dplyr::filter(gwas_olaps, gene_name %like% "FADS2", trait == "RA")

#Extract QTL
qtl_row = dplyr::filter(trqtl_min_pvalues$revisedAnnotations$AcLDL, phenotype_id == "ENSG00000134824.clique_1.upstream.ENST00000355484") %>%
  dplyr::mutate(phenotype_id = "ENSG00000134824.clique_1.upstream.ENST00000521849")
qtl_data = extractTrQtlDataFromSE(qtl_row, se_revised, vcf_file$genotypes, variant_information, phenotype_name = "tpms")
plotQtlRow(qtl_data)

#Extract all alternative transcripts in this clique
alt_promoters = dplyr::filter(gene_meta_list$revisedAnnotation, gene_id == "ENSG00000134824.upstream", transcript_id %like% "clique_1")
alt_granges = leafcutter_granges[alt_promoters$transcript_id]

#Make a coverage plot
track_data = wiggleplotrGenotypeColourGroup(str2_df, qtl_row$snp_id, vcf_file$genotypes, 1)
utr = GRanges(seqnames = "11", IRanges(start = 61846591, end =61848158 - 100), strand = "+")
lipa_coverage = wiggleplotr::plotCoverage(list(UTR = utr), track_data = track_data, 
                                          fill_palette = getGenotypePalette(), 
                                          plot_fraction = 0.2, 
                                          coverage_type = "line", 
                                          rescale_introns = FALSE, flanking_length = c(25,25))



#Make eQTL boxplot
gene_data = constructQtlPlotDataFrame("ENSG00000134824", "rs61896141", acldl_list$cqn, vcf_file$genotypes, 
                                      acldl_list$sample_metadata, acldl_list$gene_metadata) %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("Ctrl", "AcLDL"))) %>%
  dplyr::left_join(constructGenotypeText("rs61896141", variant_information), by = "genotype_value")
qtl_plot = plotQtlRow(gene_data)
ggsave("acLDL_figures/FADS2_boxplot.pdf", plot = qtl_plot, width = 4, height = 3)
