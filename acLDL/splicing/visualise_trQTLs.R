library("wiggleplotr")
library("GenomicFeatures")

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

track_data = wiggleplotrGenotypeColourGroup(str1_df, "rs12794886", vcf_file$genotypes, 1)

#Extract transcripts of the gene
selected_transcripts = dplyr::filter(gene_metadata, gene_id == "ENSG00000149541")$transcript_id
plotCoverage(exons[selected_transcripts], cdss[intersect(selected_transcripts, names(cdss))], 
             track_data, gene_metadata, 
             heights = c(2,1), fill_palette = getGenotypePalette(), coverage_type = "line", rescale_introns = FALSE)


data = constructQtlPlotDataFrame("ENST00000531383", "rs12794886", assays(se_ensembl)$tpm_ratios, vcf_file$genotypes, 
                                 tbl_df2(colData(se_ensembl)), 
                                 tbl_df2(rowData(se_ensembl)) %>% dplyr::mutate(gene_id = transcript_id)) %>%
  dplyr::left_join(constructGenotypeText("rs12794886", variant_information), by = "genotype_value")
plotQtlRow(data)


