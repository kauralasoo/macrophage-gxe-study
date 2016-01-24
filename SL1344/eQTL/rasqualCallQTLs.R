library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")

#Import all results after eigenMT correction
naive_eigenMT = eigenMTImportResults("results/SL1344/rasqual/output/naive_500kb/naive_500kb.eigenMT.txt")
IFNg_eigenMT = eigenMTImportResults("results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.eigenMT.txt")
SL1344_eigenMT = eigenMTImportResults("results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.eigenMT.txt")
IFNg_SL1344_eigenMT = eigenMTImportResults("results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.eigenMT.txt")
shared_genes = intersect(intersect(intersect(naive_eigenMT$gene_id, IFNg_eigenMT$gene_id), SL1344_eigenMT$gene_id), IFNg_eigenMT$gene_id)

#Find hits in all conditions
min_pvalue_list = list(naive = naive_eigenMT, IFNg = IFNg_eigenMT, SL1344 = SL1344_eigenMT, IFNg_SL1344 = IFNg_SL1344_eigenMT)
min_pvalue_list = lapply(min_pvalue_list, function(x, shared_genes){res = dplyr::filter(x, gene_id %in% shared_genes)}, shared_genes)
min_pvalue_hits = lapply(min_pvalue_list, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = ldply(min_pvalue_hits, .id = "condition_name")
joint_pairs = dplyr::select(min_pvalue_df, gene_id, snp_id) %>% unique() %>%
  dplyr::filter(gene_id %in% shared_genes)
saveRDS(min_pvalue_list, "results/SL1344/eQTLs/rasqual_min_pvalues.rds")

#Make database connections
naive_sql <- src_sqlite("~/projects/macrophage-gxe-study/databases/naive_500kb.sqlite3") %>% tbl("rasqual")
sl1344_sql <- src_sqlite("~/projects/macrophage-gxe-study/databases/SL1344_500kb.sqlite3") %>% tbl("rasqual")
ifng_sql <- src_sqlite("~/projects/macrophage-gxe-study/databases/IFNg_500kb.sqlite3") %>% tbl("rasqual")
ifng_sl1344_sql <- src_sqlite("~/projects/macrophage-gxe-study/databases/IFNg_SL1344_500kb.sqlite3") %>% tbl("rasqual")

#Extract p-values for all peaks that had significant QTLs
naive_pvalues = fetchMultipleGenes(joint_pairs, naive_sql)
ifng_pvalues = fetchMultipleGenes(joint_pairs, ifng_sql)
sl1344_pvalues = fetchMultipleGenes(joint_pairs, sl1344_sql)
ifng_sl1344_pvalues = fetchMultipleGenes(joint_pairs, ifng_sl1344_sql)

selected_pvalue_list = list(naive = naive_pvalues,
                            IFNg = ifng_pvalues,
                            SL1344 = sl1344_pvalues,
                            IFNg_SL1344 = ifng_sl1344_pvalues)
saveRDS(selected_pvalue_list, "results/SL1344/eQTLs/rasqual_selected_pvalues.rds")
