library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("readr")
library("devtools")
library("optparse")
load_all("../seqUtils/")


gwas_id = "IBD"
cis_window = 2e5

#Import variant information
GRCh38_variants = importVariantInformation("genotypes/acLDL/imputed_20151005/imputed.70_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/acLDL/imputed_20151005/GRCh37/imputed.70_samples.variant_information.GRCh37.txt.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name"))

#Set intput and output files
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path("databases/GWAS/summary", gwas_file_name)

#Import lead p-values for QTLs
min_pvalues = list(Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/Ctrl.permuted.txt.gz"), 
                   AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/AcLDL.permuted.txt.gz")) %>%
  purrr::map(~dplyr::select(., phenotype_id, snp_id, p_fdr))

#QTL summary stats
qtl_summary_list = list(Ctrl = "processed/acLDL/fastqtl_output/ensembl_87/sorted/Ctrl.nominal.sorted.txt.gz",
                        AcLDL = "processed/acLDL/fastqtl_output/ensembl_87/sorted/AcLDL.nominal.sorted.txt.gz")
sample_sizes = list(Ctrl = 70, AcLDL = 70)

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, 
                                       overlap_dist = 1e5, gwas_thresh = 1e-5)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()

#Test for coloc
coloc_res_list = purrr::map2(qtl_summary_list, sample_sizes, 
                             ~colocMolecularQTLsByRow(qtl_pairs[1:5,], qtl_summary_path = .x, 
                                                      gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants = GRCh37_variants, 
                                                                                      GRCh38_variants = GRCh38_variants, N_qtl = .y, cis_dist = cis_window))

#Debugging example
#colocMolecularQTLs(qtl_pairs[1,], qtl_summary_list$Ctrl, gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants, GRCh38_variants, N_qtl = 2e5)
