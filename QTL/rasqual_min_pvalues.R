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
saveRDS("results/ATAC/QTLs/rasqual_min_pvalues.txt")