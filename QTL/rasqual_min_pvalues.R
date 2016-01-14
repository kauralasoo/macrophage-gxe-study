library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")

#Extract minimal p-value for each condition
naive_min_pvalues = findMinimalSnpPvaluesSQLite("~/projects/macrophage-chromatin/naive_50kb.sqlite3")
ifng_min_pvalues = findMinimalSnpPvaluesSQLite("~/projects/macrophage-chromatin/IFNg_50kb.sqlite3")
sl1344_min_pvalues = findMinimalSnpPvaluesSQLite("~/projects/macrophage-chromatin/SL1344_50kb.sqlite3")
ifng_sl1344_min_pvalues = findMinimalSnpPvaluesSQLite("~/projects/macrophage-chromatin/IFNg_SL1344_50kb.sqlite3")

min_pvalue_list = list(naive = naive_min_pvalues,
                       IFNg = ifng_min_pvalues,
                       SL1344 = sl1344_min_pvalues,
                       IFNg_SL1344 = ifng_sl1344_min_pvalues)
saveRDS(min_pvalue_list, "results/ATAC/QTLs/rasqual_min_pvalues.rds")


#Make database connections
naive_sql <- src_sqlite("~/projects/macrophage-chromatin/naive_50kb.sqlite3") %>% tbl("rasqual")
ifng_sql <- src_sqlite("~/projects/macrophage-chromatin/IFNg_50kb.sqlite3") %>% tbl("rasqual")
sl1344_sql <- src_sqlite("~/projects/macrophage-chromatin/SL1344_50kb.sqlite3") %>% tbl("rasqual")
ifng_sl1344_sql <- src_sqlite("~/projects/macrophage-chromatin/IFNg_SL1344_50kb.sqlite3") %>% tbl("rasqual")

#Extract all SNP-gene pairs
min_pvalues_list = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalues_df = ldply(min_pvalues_hits, .id = "condition_name")
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()


#Extract p-values for all peaks that had significant QTLs
gene_ids = unique(joint_pairs$gene_id)
naive_pvalues = fetchMultipleGenes(gene_ids, naive_sql)
ifng_pvalues = fetchMultipleGenes(gene_ids, ifng_sql)
sl1344_pvalues = fetchMultipleGenes(gene_ids, sl1344_sql)
ifng_sl1344_pvalues = fetchMultipleGenes(gene_ids, ifng_sl1344_sql)
selected_pvalue_list = list(naive = naive_pvalues,
                       IFNg = ifng_pvalues,
                       SL1344 = sl1344_pvalues,
                       IFNg_SL1344 = ifng_sl1344_pvalues)
saveRDS(selected_pvalue_list, "results/ATAC/QTLs/rasqual_selected_pvalues.rds")
