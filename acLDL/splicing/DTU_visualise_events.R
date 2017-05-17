library("dplyr")
library("DRIMSeq")
library("devtools")
library("optparse")
library("wiggleplotr")
library("SummarizedExperiment")
load_all("../seqUtils/")

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

#Import revisedAnnotations Granges
revised_granges = readRDS("results/reviseAnnotations/reviseAnnotations.GRangesList.rds") %>%
  purrr::flatten()
leafcutter_granges = readRDS("results/acLDL/acLDL_leafcutter.GRangesList.rds")

#Set up sample coverage df
str1_df = wiggleplotrConstructMetadata(acldl_list$counts, 
                                       sample_meta_list$revisedAnnotation, 
                                       "/Volumes/JetDrive/bigwig/", 
                                       bigWig_suffix = ".str1.bw",
                                       condition_name_levels = c("Ctrl","AcLDL")) %>%
  dplyr::mutate(colour_group = track_id)
str2_df = wiggleplotrConstructMetadata(acldl_list$counts, 
                                       sample_meta_list$revisedAnnotation, 
                                       "/Volumes/JetDrive/bigwig/", 
                                       bigWig_suffix = ".str2.bw",
                                       condition_name_levels = c("Ctrl","AcLDL")) %>%
  dplyr::mutate(colour_group = track_id)

#Import differential events
dtu_list = readRDS("acLDL_figures/tables/DTU_genes.rds")
emsembl_metadata = colData(se_ensembl)

#Identify genes
ensembl_hits = dplyr::filter(dtu_list$ensembl_87, max_diff > 0.1, p_fdr < 0.01) %>% 
  dplyr::arrange(-max_diff) %>%
  dplyr::left_join(dplyr::select(gene_meta_list$ensembl_87, gene_id, gene_name, strand) %>% unique(), by = "gene_id")

gene_meta = dplyr::filter(gene_meta_list$ensembl_87, gene_id == "ENSG00000103569")
plotCoverage(exons[gene_meta$transcript_id], track_data = str2_df)

#Identify genes
leafcutter_hits = dplyr::filter(dtu_list$leafcutter, max_diff > 0.1, p_fdr < 0.01) %>% 
  dplyr::arrange(-max_diff) %>%
  dplyr::left_join(dplyr::select(gene_meta_list$leafcutter, gene_id, gene_name, ensembl_gene_strand) %>% unique(), by = "gene_id")

gene_meta = dplyr::filter(gene_meta_list$leafcutter, gene_id == "clu_28487")
plotCoverage(leafcutter_granges[gene_meta$transcript_id], track_data = str2_df, plot_fraction = 0.2)

#Identify genes
revised_hits = dplyr::filter(dtu_list$revisedAnnotation, max_diff > 0.1, p_fdr < 0.01) %>% 
  dplyr::arrange(-max_diff) %>%
  dplyr::left_join(dplyr::select(gene_meta_list$revisedAnnotation, gene_id, gene_name, strand) %>% unique(), by = "gene_id")

gene_meta = dplyr::filter(gene_meta_list$revisedAnnotation, gene_id == "ENSG00000166523.downstream")
plotCoverage(revised_granges[gene_meta$transcript_id], track_data = str1_df)




