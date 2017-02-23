library("wiggleplotr")
library("GenomicFeatures")

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
variant_information = importVariantInformation("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

#Import tranqscript expression data
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
gene_metadata = rowData(se_ensembl) %>% tbl_df2()
sample_metadata = colData(se_ensembl) %>% tbl_df2()

#Import transcript annotations
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

#Set up sample coverage df
str1_df = wiggleplotrConstructMetadata(acldl_list$counts, 
                                       sample_metadata, 
                                       "/Volumes/JetDrive/bigwig/", 
                                       bigWig_suffix = ".str1.bw",
                                       condition_name_levels = c("Ctrl","AcLDL"))
str2_df = wiggleplotrConstructMetadata(acldl_list$counts, 
                                       sample_metadata, 
                                       "/Volumes/JetDrive/bigwig/", 
                                       bigWig_suffix = ".str2.bw",
                                       condition_name_levels = c("Ctrl","AcLDL"))

#Import interaction results
interaction_list = readRDS("results/acLDL/trQTLs/trQTL_interaction_results.rds")

#Import Diff p-values
diff_pvals = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/Diff.permuted.txt.gz") %>%
  dplyr::filter(p_fdr < 0.1) %>%
  dplyr::transmute(gene_id = group_id, transcript_id = phenotype_id, snp_id, strand)


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
plots_df = purrr::by_row(ensembl_interactions[1:20,], 
                         ~makeQTLCoveragePlot(., str1_df, str2_df, vcf_file$genotypes, gene_metadata, exons, cdss, 
                                              heights = c(2,1), coverage_type = "line", rescale_introns = TRUE),
                         .to = "plots")
#Save coverage plots to disk
plot_list = plots_df$plots
names(plot_list) = plots_df$transcript_id
savePlotList(plot_list, "results/acLDL/trQTLs/ensembl_plots/")



