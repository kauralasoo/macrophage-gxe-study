library("devtools")
library("plyr")
library("SummarizedExperiment")
library("dplyr")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")

#Import data
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")

#### eigenMT ####
#Export genotype data
chromosome_list = scan("macrophage-gxe-study/data/sample_lists/chromosome_list.txt", what = "char")
eigenMTExportGenotypesByChr(chromosome_list, "processed/Fairfax/geno_by_chr/",
                            "processed/Fairfax/eigenMT/input/", "")

#Export gene metadata
gene_meta = rowData(se_fairfax) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id) %>%
  dplyr::rename(start = gene_start, end = gene_end)
eigenMTExportGeneMetadata(gene_meta, "processed/Fairfax/eigenMT/input/")


#Export eQTLs in eigenMT format
qtl_min_pvalues = readRDS("results/Fairfax/fairfax_qtl_min_pvalues.rds")

#Make QTL dfs
qtl_df_list = purrr::map(qtl_min_pvalues$full, ~dplyr::transmute_(., 'SNP' = 'snp_id', 'gene' = 'phenotype_id', 'beta' = 'slope', 't-stat'= 1, 'p-value' = 'p_nominal', 'FDR' = 'p_fdr'))

#Make file names list
file_names = paste0("processed/Fairfax/eigenMT/input/", names(qtl_df),".input.txt")
file_names_list = as.list(file_names)
names(file_names_list) = names(qtl_df)

#Export QTLs
purrr::map2(qtl_df_list, file_names_list, ~write.table(.x,.y, quote = FALSE, row.names = FALSE, sep = "\t"))
