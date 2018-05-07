library("dplyr")
library("tidyr")
library("purrr")
library("UpSetR")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("../wiggleplotr/")

#Import SummarizedExperiment
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")
expression_mat = assays(se_fairfax)$exprs
sample_metadata = colData(se_fairfax) %>% tbl_df2()
gene_metadata = rowData(se_fairfax) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id)

#Import genotypes
vcf_file = seqUtils::gdsToMatrix("processed/Fairfax/geno_by_chr/20.gds")
variant_information = importVariantInformation("processed/Fairfax/merged_genotypes/fairfax_genotypes.variant_information.txt.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name","type"))

#Make a QTL plot
plot_data = constructQtlPlotDataFrame("ILMN_1779257", "rs4239702", 
                                       expression_mat, 
                                       vcf_file$genotypes, 
                                       sample_metadata, 
                                       gene_metadata) %>%
  dplyr::left_join(constructGenotypeText("rs6074021", variant_information), by = "genotype_value")

plotQtlCol(plot_data, scales = "free_y")


#Import QTL and GWAS summary stats and convert them to the same GRCh38 coordinate space
qtl_df = data_frame(phenotype_id = "ILMN_1779257", snp_id = "rs4239702", trait = "RA")
qtl_paths = list(CD14 = "~/databases/Fairfax/full/CD14.nominal.sorted.txt.gz",
                 IFN = "~/databases/Fairfax/full/IFN.nominal.sorted.txt.gz",
                 LPS2 = "~/databases/Fairfax/full/LPS2.nominal.sorted.txt.gz",
                 LPS24 = "~/databases/Fairfax/full/LPS24.nominal.sorted.txt.gz")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, gwas_dir = "~/datasets/Inflammatory_GWAS/", qtl_paths = qtl_paths, 
                                         GRCh37_variants = variant_information, GRCh38_variants = variant_information, cis_dist = 2e5, type = "QTLTools") %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::rename(track_id = condition_name) %>%
  dplyr::mutate(track_id = factor(track_id, levels = c("RA","CD14","IFN","LPS2","LPS24")))


region_coords = c(min(qtl_summary$pos), max(qtl_summary$pos))
wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id %in% c("RA","IFN")), ,region_coords = region_coords, color_R2 = TRUE)
