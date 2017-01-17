library("devtools")
library("plyr")
library("dplyr")
library("ggplot2")
library("purrr")
library("GenomicFeatures")
load_all("../seqUtils/")
load_all("../wiggleplotr")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#ATAC-seq
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, 
                                                  levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")
variant_information = importVariantInformation("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

#Import old and new variant coordinates
GRCh38_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name"))

#Import transcript annotations and metadata
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_79/TranscriptDb_GRCh38_79.db")
tx_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_79/Homo_sapiens.GRCh38.79.transcript_data.rds") %>%
  dplyr::rename(transcript_id = ensembl_transcript_id,
                gene_id = ensembl_gene_id,
                gene_name = external_gene_name)
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

#Make a boxplot of the eQTL data
TRAF1_data = constructQtlPlotDataFrame("ENSG00000056558", "rs10985070", 
                                       combined_expression_data$cqn, 
                                       vcf_file$genotypes, 
                                       combined_expression_data$sample_metadata, 
                                       combined_expression_data$gene_metadata) %>%
  dplyr::left_join(constructGenotypeText("rs10985070", variant_information), by = "genotype_value") %>%
  dplyr::filter(condition_name %in% c("naive", "IFNg_SL1344")) %>%
  dplyr::left_join(figureNames(), by = "condition_name") %>%
  dplyr::mutate(condition_name = figure_name)
TRAF1_col_plot = plotQtlCol(TRAF1_data, scales = "free_y") + 
  ggplot2::scale_color_manual(values = conditionPalette()[c(1,4)], guide=FALSE)
ggsave("figures/main_figures/TRAF1_boxplot_column.pdf", TRAF1_col_plot, width = 2, height = 3.5)


#Plot gene structure
region_coords = c(120850000, 120975000)
both_tx = dplyr::filter(tx_metadata, gene_name %in% c("TRAF1"), 
                        transcript_biotype == "protein_coding", transcript_status == "KNOWN")
tx_plot = wiggleplotr::plotTranscripts(exons[both_tx$transcript_id[1]], cdss[both_tx$transcript_id[1]], 
                                       tx_metadata, rescale_introns = FALSE, region_coords = region_coords) +
  scale_x_continuous(limits = region_coords, expand = c(0,0), breaks=pretty_breaks(n=4))


#Import QTL and GWAS summary stats and convert them to the same GRCh38 coordinate space
qtl_df = data_frame(gene_id = "ENSG00000056558", snp_id = "rs10985070", trait = "RA")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, qtl_paths = qtlResults()$rna_fastqtl, 
                           GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants, cis_dist = 2e5) %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::filter(condition_name %in% c("RA","naive","IFNg_SL1344")) %>%
  dplyr::mutate(condition_name = as.vector(condition_name)) %>% 
  dplyr::mutate(condition_name = dplyr::case_when(.$condition_name == "naive" ~ "N", 
                                                  .$condition_name == "IFNg_SL1344" ~ "I+S", 
                                                  TRUE ~ .$condition_name)) %>%
  dplyr::rename(track_id = condition_name) %>%
  dplyr::mutate(track_id = factor(track_id, levels = c("RA","N","I+S")))

#Make a manhattan plot
qtl_manhattan = makeManhattanPlot(qtl_summary, region_coords, color_R2 = TRUE)

#Make a joint plot
joint_plot = cowplot::plot_grid(qtl_manhattan, tx_plot, 
                                align = "v", ncol = 1, rel_heights = c(4,1.5))
ggsave("figures/main_figures/TRAF1_manhattan_plot.pdf", joint_plot, width = 4, height = 4)


#See if there are any associated ATAC peaks (apparently not)
a = tabixFetchSNPsQuick("rs4511811", qtlResults()$atac_rasqual$naive, vcf_file$snpspos)
a = tabixFetchSNPsQuick("rs10985070", qtlResults()$atac_rasqual$IFNg, vcf_file$snpspos)

