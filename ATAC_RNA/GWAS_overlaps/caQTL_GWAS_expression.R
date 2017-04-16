library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import eQTL p-values
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
eigenMT_tests = purrr::map_df(rasqual_min_pvalues, ~dplyr::select(., gene_id, n_tests)) %>%
  unique()

fastqtl_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")
fastqtl_lead_pvalues_df = purrr::map_df(fastqtl_min_pvalues, identity, .id = "condition_name") %>%
  dplyr::select(condition_name, snp_id, gene_id) %>%
  dplyr::rename(eQTL_snp_id = snp_id)
fastqtl_n_snps = purrr::map_df(fastqtl_min_pvalues, ~dplyr::select(., gene_id, n_cis_snps)) %>%
  unique()

#Lead SNP df
#lead_pvalues_df = purrr::map_df(rasqual_min_pvalues, identity, .id = "condition_name") %>%
#  dplyr::select(condition_name, snp_id, gene_id) %>%
#  dplyr::rename(eQTL_snp_id = snp_id)

#Import colocalised caQTLs
caqtl_200kb_filtered_hits = readRDS("results/SL1344/coloc/caQTL_coloc_200kb_hits.rds")

#Identify unique GWAS hits
unique_hits = dplyr::select(caqtl_200kb_filtered_hits,phenotype_id, snp_id) %>% unique()
unique_hits_traits = dplyr::select(caqtl_200kb_filtered_hits,phenotype_id, snp_id, trait) %>% unique()

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#for each SNP, fetch p-values for all tested genes
selected_snps = dplyr::filter(vcf_file$snpspos, snpid %in% unique_hits$snp_id) %>%
  dplyr::transmute(snp_id = snpid, seqnames = chr, start = pos, end = pos, strand = "*") %>%
  dataFrameToGRanges()

snp_pvalues = purrr::map_df(qtlResults()$rna_fastqtl, ~fastqtlTabixFetchSNPs(selected_snps, .), .id = "condition_name") %>% 
  dplyr::rename(gene_id = phenotype_id)

filtered_overlaps = dplyr::left_join(snp_pvalues, fastqtl_n_snps, by = "gene_id") %>%
  dplyr::mutate(p_eigen = pmin(p_nominal*n_cis_snps,1)) %>%
  dplyr::filter(!is.na(p_eigen)) %>%
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(gene_id, p_eigen) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(p_eigen < 0.1) %>%
  dplyr::select(gene_id, snp_id, condition_name, p_nominal, p_eigen) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") %>%
  dplyr::left_join(unique_hits_traits, by = "snp_id") %>%
  dplyr::left_join(fastqtl_lead_pvalues_df, by = c("gene_id","condition_name")) %>%
  dplyr::group_by(gene_id, snp_id) %>%
  dplyr::mutate(R2 = calculatePairR2(snp_id, eQTL_snp_id, vcf_file$genotypes))
write.table(filtered_overlaps, "results/ATAC_RNA_overlaps/caQTL_eQTL_rescued_pairs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Import old and new variant coordinates
GRCh38_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")



#Import eQTL summary stats for the ICOSLG gene



plotColocs <- function(gene_id, peak_id, snp, trait_id, rna_paths = qtlResults()$rna_rasqual, rna_type = "RASQUAL"){
   qtl_df = data_frame(phenotype_id = peak_id, snp_id = snp, trait = trait_id)
  qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, qtl_paths = qtlResults()$atac_fastqtl, 
                                           GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants, 
                                           cis_dist = 2e5, gwas_dir = "~/datasets/Inflammatory_GWAS/") %>%
    arrange(condition_name, p_nominal) %>% 
    addR2FromLead(vcf_file$genotypes) %>%
    dplyr::mutate(track_id = condition_name)
  region_coords = c(min(qtl_summary$pos), max(qtl_summary$pos))
  caqtl_manhattan_plot = makeManhattanPlot(qtl_summary, region_coords, color_R2 = TRUE, data_track = TRUE)
  
  
  #Import eQTL summary stats for the ICOSLG gene
  qtl_df = data_frame(phenotype_id = gene_id, snp_id = snp, trait = trait_id)
  eqtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, qtl_paths = rna_paths, 
                                            GRCh37_variants = GRCh37_variants, GRCh38_variants = GRCh38_variants,
                                            cis_dist = 2e5, gwas_dir = "~/datasets/Inflammatory_GWAS/", type = rna_type) %>%
    arrange(condition_name, p_nominal) %>% 
    addR2FromLead(vcf_file$genotypes) %>%
    dplyr::mutate(track_id = condition_name)
  region_coords = c(min(qtl_summary$pos), max(qtl_summary$pos))
  eqtl_manhattan_plot = makeManhattanPlot(eqtl_summary, region_coords, color_R2 = TRUE, data_track = TRUE)
  
  return(list(caqtl_manhattan_plot,eqtl_manhattan_plot))
}

all_plots = purrr::by_row(filtered_overlaps, ~plotColocs(.$gene_id, .$phenotype_id, .$snp_id, .$trait))
all_rasqual = purrr::by_row(filtered_overlaps, ~plotColocs(.$gene_id, .$phenotype_id, .$snp_id, .$trait))



#Explore the CTSB locus in detail
gene_data = constructQtlPlotDataFrame("ENSG00000164733", "rs11997338", combined_expression_data$cqn, vcf_file$genotypes, 
                                      combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg_SL1344")) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::left_join(constructGenotypeText("rs11997338", variant_information), by = "genotype_value")
gene_plot = plotQTLCompact(gene_data) + ggplot2::scale_color_manual(values = conditionPalette()[c(1,4)], guide=FALSE)
ggsave("figures/main_figures/CTSC_expression_boxplot.pdf", plot = gene_plot, width = 2, height = 2.5)

peak_data = constructQtlPlotDataFrame("ATAC_peak_260182", "rs11997338", atac_list$cqn, vcf_file$genotypes, 
                                      atac_list$sample_metadata, atac_list$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg_SL1344")) %>%
  dplyr::left_join(figureNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::left_join(constructGenotypeText("rs11997338", variant_information), by = "genotype_value")
peak_plot = plotQTLCompact(peak_data) + ggplot2::scale_color_manual(values = conditionPalette()[c(1,4)], guide=FALSE) + 
  ylab("Chromatin accessibility")


