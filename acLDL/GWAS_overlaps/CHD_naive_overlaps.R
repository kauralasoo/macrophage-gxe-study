library("dplyr")
library("tidyr")
library("SummarizedExperiment")

#Import variant information
GRCh38_variants = importVariantInformation("genotypes/acLDL/imputed_20151005/imputed.70_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/acLDL/imputed_20151005/GRCh37/imputed.70_samples.variant_information.GRCh37.txt.gz")

#Import the VCF file
vcf_file = readRDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.rds")

#Import expression data
combined_expression_data = readRDS("results/acLDL/acLDL_combined_expression_data.rds")
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name) %>% 
  dplyr::rename(phenotype_id = gene_id)

#Load eQTL p-values from disk
rasqual_min_pvalues = readRDS("results/acLDL/eQTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Import CHD hits
chd_hits = read.table("macrophage-gxe-study/data/gwas_catalog/UKBB_Exome_CHD_Loci.txt", 
                      header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::tbl_df() %>%
  tidyr::separate("Chr.Pos", into = c("chr","pos"), sep = ":") %>%
  dplyr::mutate(snp_id = Markername) %>%
  dplyr::mutate(pos = as.integer(pos)) %>%
  summaryReplaceSnpId(GRCh37_variants) %>% 
  summaryReplaceCoordinates(GRCh38_variants)

#eQTL overlaps
chd_eQTL_overlaps = findGWASOverlaps(joint_pairs, chd_hits, vcf_file, max_distance = 5e5, min_r2 = 0.8) %>% 
  dplyr::left_join(gene_name_map, by = c("gene_id" = "phenotype_id")) %>%
  dplyr::select(gene_id, Markername, R2, gene_name) %>% dplyr::group_by(Markername, gene_name) %>%
  arrange(Markername, gene_name, -R2) %>%
  dplyr::filter(row_number() == 1)  %>%
  dplyr::mutate(phenotype = "eQTL_rasqual")

#Do a quick scan for trQTLs as well
phenotype_list = list(
  #Ensembl 87 trQTLs
  ensembl_87 = list(
    min_pvalues = list(Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/Ctrl.permuted.txt.gz"), 
                       AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/AcLDL.permuted.txt.gz")) %>%
      purrr::map(~dplyr::select(., phenotype_id, snp_id, p_fdr)),
    qtl_summary_list = list(Ctrl = "processed/acLDL/fastqtl_output/ensembl_87/sorted/Ctrl.nominal.sorted.txt.gz",
                            AcLDL = "processed/acLDL/fastqtl_output/ensembl_87/sorted/AcLDL.nominal.sorted.txt.gz"),
    sample_sizes = list(Ctrl = 70, AcLDL = 70)
  ),
  reviseAnnotations = list(
    min_pvalues = list(Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/reviseAnnotations/Ctrl.permuted.txt.gz"), 
                       AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/reviseAnnotations/AcLDL.permuted.txt.gz")) %>%
      purrr::map(~dplyr::select(., phenotype_id, snp_id, p_fdr)),
    qtl_summary_list = list(Ctrl = "processed/acLDL/fastqtl_output/reviseAnnotations/sorted/Ctrl.nominal.sorted.txt.gz",
                            AcLDL = "processed/acLDL/fastqtl_output/reviseAnnotations/sorted/AcLDL.nominal.sorted.txt.gz"),
    sample_sizes = list(Ctrl = 70, AcLDL = 70)
  ),
  leafcutter = list(
    min_pvalues = list(Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/leafcutter/Ctrl.permuted.txt.gz"), 
                       AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/leafcutter/AcLDL.permuted.txt.gz")) %>%
      purrr::map(~dplyr::select(., phenotype_id, snp_id, p_fdr)),
    qtl_summary_list = list(Ctrl = "processed/acLDL/fastqtl_output/leafcutter/sorted/Ctrl.nominal.sorted.txt.gz",
                            AcLDL = "processed/acLDL/fastqtl_output/leafcutter/sorted/AcLDL.nominal.sorted.txt.gz"),
    sample_sizes = list(Ctrl = 70, AcLDL = 70)
  )
)

#Import SummarizedExperiments
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
se_reviseAnnotations = readRDS("results/acLDL/acLDL_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/acLDL/acLDL_leafcutter_counts.rds")

#Gene names
ensembl_name_map = dplyr::select(tbl_df2(rowData(se_ensembl)), transcript_id, gene_name) %>% dplyr::rename(phenotype_id = transcript_id)
revised_name_map = dplyr::select(tbl_df2(rowData(se_reviseAnnotations)), transcript_id, gene_name) %>% dplyr::rename(phenotype_id = transcript_id)
leafcutter_name_map = dplyr::select(tbl_df2(rowData(se_leafcutter)), transcript_id, ensembl_gene_id) %>% 
  dplyr::left_join(dplyr::transmute(tbl_df2(rowData(se_ensembl)), ensembl_gene_id = gene_id, gene_name), by = "ensembl_gene_id") %>%
  dplyr::rename(phenotype_id = transcript_id) %>%
  unique()

ensembl_min_p = purrr::map_df(phenotype_list$ensembl_87$min_pvalues, identity, .id = "condition_name") %>%
  dplyr::arrange(phenotype_id, p_fdr) %>%
  dplyr::filter(p_fdr < 0.1) %>%
  dplyr::select(phenotype_id, snp_id) %>% 
  unique() %>%
  dplyr::rename(gene_id = phenotype_id)

chd_trQTL_overlaps = findGWASOverlaps(ensembl_min_p, chd_hits, vcf_file, max_distance = 5e5, min_r2 = 0.8) %>%
  dplyr::group_by(gene_id, Markername) %>%
  dplyr::arrange(gene_id, Markername, -R2) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::select(gene_id, Markername, R2) %>%
  ungroup() %>%
  dplyr::left_join(ensembl_name_map, by = c("gene_id" = "phenotype_id"))  %>%
  dplyr::mutate(phenotype = "trQTL_ensembl")
  

revised_min_p = purrr::map_df(phenotype_list$reviseAnnotations$min_pvalues, identity, .id = "condition_name") %>%
  dplyr::arrange(phenotype_id, p_fdr) %>%
  dplyr::filter(p_fdr < 0.1) %>%
  dplyr::select(phenotype_id, snp_id) %>% 
  unique() %>%
  dplyr::rename(gene_id = phenotype_id)

chd_revised_overlaps = findGWASOverlaps(revised_min_p, chd_hits, vcf_file, max_distance = 5e5, min_r2 = 0.8) %>%
  dplyr::group_by(gene_id, Markername) %>%
  dplyr::arrange(gene_id, Markername, -R2) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::select(gene_id, Markername, R2) %>%
  ungroup() %>%
  dplyr::left_join(revised_name_map, by = c("gene_id" = "phenotype_id")) %>%
  dplyr::mutate(phenotype = "trQTL_revised")

all_overlaps = dplyr::bind_rows(chd_eQTL_overlaps, chd_trQTL_overlaps, chd_revised_overlaps)
write.table(all_overlaps,"results/acLDL/naive_CHD_overlaps.txt", sep = "\t", quote = FALSE, row.names = FALSE)
