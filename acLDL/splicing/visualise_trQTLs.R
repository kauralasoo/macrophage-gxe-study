library("devtools")
library("rtracklayer")
library("wiggleplotr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../reviseAnnotations/")

#Functions
makeQTLCoveragePlot <- function(qtl_df, str1_df, str2_df, genotypes, gene_metadata, exons, cdss, ...){
  
  #Construct track data
  if(qtl_df$strand == "-"){
    bigwig_meta = str1_df
  } else if (qtl_df$strand == -1){
    bigwig_meta = str1_df
  }
  else{
    bigwig_meta = str2_df
  }
  track_data = wiggleplotrGenotypeColourGroup(bigwig_meta, qtl_df$snp_id, genotypes, 1)
  
  #Select transcripts
  selected_transcripts = dplyr::filter(gene_metadata, gene_id == qtl_df$gene_id)$transcript_id
  
  #Make a coverage plot
  plot = plotCoverage(exons[selected_transcripts], 
                      cdss[intersect(selected_transcripts, names(cdss))], 
                      track_data, gene_metadata, 
                      fill_palette = getGenotypePalette(), ...)
  return(plot)
}

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

#Import revisedAnnotations Granges
revised_granges = readRDS("results/reviseAnnotations/reviseAnnotations.GRangesList.rds")
granges = purrr::flatten(revised_granges)

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

#Import interaction results
interaction_list = readRDS("results/acLDL/trQTLs/trQTL_interaction_results.rds")


#CD33 example
cd33_hit = dplyr::filter(trqtl_min_pvalues$revisedAnnotations$Ctrl, group_id == "ENSG00000105383.contained")
tx_names = c(cd33_hit$phenotype_id, cd33_hit$other_phenotype_id)
exons = granges[tx_names]
track_data = wiggleplotrGenotypeColourGroup(str2_df, cd33_hit$snp_id, vcf_file$genotypes, 1)
wiggleplotr::plotCoverage(exons, track_data = track_data, fill_palette = getGenotypePalette(), 
                          plot_fraction = 0.2, coverage_type = "line", rescale_introns = TRUE)



#Make coverage plots for all genes
plots_df = purrr::by_row(diff_pvals, 
              ~makeQTLCoveragePlot(., str1_df, str2_df, vcf_file$genotypes, gene_metadata, exons, cdss, 
                                              heights = c(2,1), coverage_type = "line", rescale_introns = TRUE),
              .to = "plots")

#Save coverage plots to disk
plot_list = plots_df$plots
names(plot_list) = plots_df$transcript_id
savePlotList(plot_list, "results/acLDL/trQTLs/diff_plots/")


data = constructQtlPlotDataFrame("ENST00000397147", "rs2072711", assays(se_ensembl)$tpm_ratios, vcf_file$genotypes, 
                                 tbl_df2(colData(se_ensembl)), 
                                 tbl_df2(rowData(se_ensembl)) %>% dplyr::mutate(gene_id = transcript_id)) %>%
  dplyr::left_join(constructGenotypeText("rs2072711", variant_information), by = "genotype_value")
plotQtlRow(data)



data = constructQtlPlotDataFrame("ENST00000557352", "rs11394080", ensembl_abundances, vcf_file$genotypes, 
                                 tbl_df2(colData(se_ensembl)), 
                                 tbl_df2(rowData(se_ensembl)) %>% dplyr::mutate(gene_id = transcript_id)) %>%
  dplyr::left_join(constructGenotypeText("rs111343454", variant_information), by = "genotype_value")
plotQtlRow(data)


data = constructQtlPlotDataFrame("ENST00000557352", "rs11394080", assays(se_ensembl)$tpm_ratios, vcf_file$genotypes, 
                                 tbl_df2(colData(se_ensembl)), 
                                 tbl_df2(rowData(se_ensembl)) %>% dplyr::mutate(gene_id = transcript_id)) %>%
  dplyr::left_join(constructGenotypeText("rs11879855", variant_information), by = "genotype_value")
plotQtlRow(data)

ensembl_interactions = dplyr::rename(interaction_list$Ensembl, transcript_id = gene_id) %>% 
  dplyr::left_join(dplyr::select(gene_metadata, transcript_id, gene_id, gene_name, strand), by = "transcript_id") %>%
  dplyr::filter(p_fdr < 0.01)

#Make coverage plots for all genes
plots_df = purrr::by_row(ensembl_interactions, 
                         ~makeQTLCoveragePlot(., str1_df, str2_df, vcf_file$genotypes, gene_metadata, exons, cdss, 
                                              heights = c(2,1), coverage_type = "line", rescale_introns = TRUE),
                         .to = "plots")
#Save coverage plots to disk
plot_list = plots_df$plots
names(plot_list) = plots_df$transcript_id
savePlotList(plot_list, "results/acLDL/trQTLs/ensembl_plots/")


#APOBR example
data = constructQtlPlotDataFrame("ENST00000564831", "rs149271", assays(se_ensembl)$tpm_ratios, vcf_file$genotypes, 
                                 tbl_df2(colData(se_ensembl)), 
                                 tbl_df2(rowData(se_ensembl)) %>% dplyr::mutate(gene_id = transcript_id)) %>%
  dplyr::left_join(constructGenotypeText("rs149271", variant_information), by = "genotype_value")
plotQtlRow(data)

gene_info = dplyr::filter(ctrl_pvals, transcript_id == "ENST00000564831")
#Make coverage plots for all genes
plots_df = purrr::by_row(gene_info, 
                         ~makeQTLCoveragePlot(., str1_df, str2_df, vcf_file$genotypes, gene_metadata, exons, cdss, 
                                              heights = c(2,1), coverage_type = "line", rescale_introns = FALSE),
                         .to = "plots")


#PARP12
data = constructQtlPlotDataFrame("ENST00000491598", "rs7805521", assays(se_ensembl)$tpm_ratios, vcf_file$genotypes, 
                                 tbl_df2(colData(se_ensembl)), 
                                 tbl_df2(rowData(se_ensembl)) %>% dplyr::mutate(gene_id = transcript_id)) %>%
  dplyr::left_join(constructGenotypeText("rs149271", variant_information), by = "genotype_value")
plotQtlRow(data)

data = constructQtlPlotDataFrame("ENSG00000059378.clique_1.contained.ENST00000263549", "rs117102789", assays(se_reviseAnnotations)$tpm_ratios, vcf_file$genotypes, 
                                 tbl_df2(colData(se_reviseAnnotations)), 
                                 tbl_df2(rowData(se_reviseAnnotations)) %>% dplyr::mutate(gene_id = transcript_id)) %>%
  dplyr::left_join(constructGenotypeText("rs149271", variant_information), by = "genotype_value")
plotQtlRow(data)

gene_info = dplyr::filter(ctrl_pvals, gene_id == "ENSG00000105383")
#Make coverage plots for all genes
plots_df = purrr::by_row(gene_info, 
                         ~makeQTLCoveragePlot(., str1_df, str2_df, vcf_file$genotypes, gene_metadata, exons, cdss, 
                                              heights = c(2,1), coverage_type = "line", rescale_introns = TRUE),
                         .to = "plots")
plots_df$plots

