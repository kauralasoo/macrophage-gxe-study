naive_pvalues = read.table("results/SL1344/fastqtl/output/naive_filtered.txt", stringsAsFactors = FALSE)
colnames(naive_pvalues) = c("gene_id", "CHR_ID","CHR_POS", "snp_id", "distance", "pvalue", "effect_size")
naive_pvalues = dplyr::filter(naive_pvalues, CHR_ID != "X") %>% dplyr::mutate(CHR_ID = as.numeric(CHR_ID))

gwas_catalog = readr::read_delim("../../annotations/gwas/gwas_catalog_v1.0-downloaded_2015-11-18.tsv", delim = "\t")

b = dplyr::semi_join(gwas_catalog, naive_pvalues, by = c("CHR_ID", "CHR_POS"))

load("../../annotations/gwas/snp.info.Rbin")
