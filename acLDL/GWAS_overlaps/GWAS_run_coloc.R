library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("readr")
library("devtools")
library("optparse")
load_all("../seqUtils/")


gwas_id = "IBD"

#Import variant information
GRCh38_variants = importVariantInformation("genotypes/acLDL/imputed_20151005/imputed.70_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/acLDL/imputed_20151005/GRCh37/imputed.70_samples.variant_information.GRCh37.txt.gz")


#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name"))

#Set intput and output files
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path("databases/GWAS/summary", gwas_file_name)


gwas_ranges = constructVariantRanges(qtl_df, GRCh37_variants, cis_dist = cis_dist)

min_pvalues = list(Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/Ctrl.permuted.txt.gz") %>%
                     dplyr::transmute(gene_id = phenotype_id, snp_id, p_fdr))

lead_ranges = constructVariantRanges(dplyr::select(min_pvalues[1:2,], snp_id, phenotype_id), GRCh38_variants, 2e5)
pvals = qtltoolsTabixFetchPhenotypes(lead_ranges, "processed/acLDL/fastqtl_output/ensembl_87/sorted/Ctrl.nominal.sorted.txt.gz")


qtl_df_list = prefilterColocCandidates(min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, 
                                       overlap_dist = 1e5, gwas_thresh = 1e-5)


