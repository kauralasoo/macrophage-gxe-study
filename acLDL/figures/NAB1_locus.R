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

#identify QTLs
ensembl_qtls = dplyr::filter(trqtl_min_pvalues$ensembl_87$Ctrl, group_id == "ENSG00000138386")
revised_qtls = dplyr::filter(trqtl_min_pvalues$revisedAnnotations$Ctrl, group_id %like% "ENSG00000138386")
leafcutter_qtls = dplyr::filter(trqtl_min_pvalues$leafcutter$Ctrl, group_id == "clu_8618")

#Extract QTL
downstream_qtl = revised_qtls[2,]
downstream_tx_names = dplyr::filter(gene_meta_list$revisedAnnotation, gene_id == "ENSG00000138386.downstream")$transcript_id
downstream_tx = revised_granges[downstream_tx_names]

#Extract annotatated transcripts
lipa_tx = exons["ENST00000336233"] %>% removeMetadata()
lipa_cdss = cdss["ENST00000336233"] %>% removeMetadata()

#Make a coverage plot
track_data = wiggleplotrGenotypeColourGroup(str2_df, downstream_qtl$snp_id, vcf_file$genotypes, 1)
lipa_coverage = wiggleplotr::plotCoverage(downstream_tx, track_data = track_data, 
                                          fill_palette = getGenotypePalette(), 
                                          plot_fraction = 0.2, 
                                          coverage_type = "line", 
                                          rescale_introns = TRUE)
ggsave("acLDL_figures/LIPA_leafcutter_coverage.pdf", lipa_coverage, width = 10, height = 6)


#Extract QTL
downstream_qtl = ensembl_qtls
downstream_tx_names = dplyr::filter(gene_meta_list$ensembl_87, gene_id == "ENSG00000138386")$transcript_id
downstream_tx = exons[downstream_tx_names]

#Extract annotatated transcripts
lipa_tx = exons["ENST00000336233"] %>% removeMetadata()
lipa_cdss = cdss["ENST00000336233"] %>% removeMetadata()

#Make a coverage plot
track_data = wiggleplotrGenotypeColourGroup(str2_df, downstream_qtl$snp_id, vcf_file$genotypes, 1)
lipa_coverage = wiggleplotr::plotCoverage(downstream_tx, track_data = track_data, 
                                          fill_palette = getGenotypePalette(), 
                                          plot_fraction = 0.2, 
                                          coverage_type = "line", 
                                          rescale_introns = TRUE)
ggsave("acLDL_figures/LIPA_leafcutter_coverage.pdf", lipa_coverage, width = 10, height = 6)

