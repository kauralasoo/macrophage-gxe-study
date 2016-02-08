library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("Rsamtools")

#Import data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")

#Import SNP coordinates
snp_coords = read.table("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.snp_coords.txt", stringsAsFactors = FALSE)
colnames(snp_coords) = c("chr", "pos", "snp_id")

#Extract minimal p-value for each condition
naive_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output/naive_100kb/naive_50kb.eigenMT.txt")
ifng_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output/IFNg_100kb/IFNg_50kb.eigenMT.txt")
sl1344_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output/SL1344_100kb/SL1344_50kb.eigenMT.txt")
ifng_sl1344_eigen_pvalue = eigenMTImportResults("results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_50kb.eigenMT.txt")

min_pvalue_list = list(naive = naive_eigen_pvalue,
                       IFNg = ifng_eigen_pvalue,
                       SL1344 = sl1344_eigen_pvalue,
                       IFNg_SL1344 = ifng_sl1344_eigen_pvalue)
saveRDS(min_pvalue_list, "results/ATAC/QTLs/rasqual_min_pvalues.rds")


#Make database connections
naive_sql <- src_sqlite("~/projects/macrophage-chromatin/naive_100kb.sqlite3") %>% tbl("rasqual")
ifng_sql <- src_sqlite("~/projects/macrophage-chromatin/IFNg_100kb.sqlite3") %>% tbl("rasqual")
sl1344_sql <- src_sqlite("~/projects/macrophage-chromatin/SL1344_100kb.sqlite3") %>% tbl("rasqual")
ifng_sl1344_sql <- src_sqlite("~/projects/macrophage-chromatin/IFNg_SL1344_100kb.sqlite3") %>% tbl("rasqual")

naive_sql <- src_sqlite("/Volumes/JetDrive/ATAC_sqlite/naive_100kb.sqlite3") %>% tbl("rasqual")
IFNg_sql <- src_sqlite("/Volumes/JetDrive/ATAC_sqlite/IFNg_100kb.sqlite3") %>% tbl("rasqual")


#Extract all SNP-gene pairs
min_pvalues_list = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalues_df = ldply(min_pvalues_hits, .id = "condition_name")
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Extract p-values for all peaks that had significant QTLs
naive_pvalues = fetchMultipleGeneSNPPairs(joint_pairs, naive_sql)
naive_pvalues = fetchMultipleGenes(joint_pairs[1,], naive_sql)

naive_pvalues = fetchMultipleGeneSNPPairs(joint_pairs[1,], IFNg_sql)
naive_pvalues = fetchMultipleGenes(joint_pairs[1,], IFNg_sql)

ifng_pvalues = fetchMultipleGenes(joint_pairs, ifng_sql)
sl1344_pvalues = fetchMultipleGenes(gene_ids, sl1344_sql)
ifng_sl1344_pvalues = fetchMultipleGenes(gene_ids, ifng_sl1344_sql)
selected_pvalue_list = list(naive = naive_pvalues,
                       IFNg = ifng_pvalues,
                       SL1344 = sl1344_pvalues,
                       IFNg_SL1344 = ifng_sl1344_pvalues)
saveRDS(selected_pvalue_list, "results/ATAC/QTLs/rasqual_selected_pvalues.rds")

#Use Rsamtools to extract 
granges = constructGeneRanges(top_hits, atac_list$gene_metadata, cis_window = 100000)
pvalues = tabixFetchGenes(granges, tabix_file)

#Construct credible sets for associations
naive_credible_set = lapply(pvalues, function(x, n, thresh){ addAssociationPosterior(x, n) %>% constructCredibleSet(thresh)}, 27, 0.99)

#Find credible sets for these peaks
ifng_granges = constructGeneRanges(ifng_effect_sizes, atac_list$gene_metadata, cis_window = 100000)
ifng_pvalues = tabixFetchGenes(ifng_granges, naive_tabix)
ifng_credible_set = lapply(ifng_pvalues, function(x, n, thresh){ addAssociationPosterior(x, n) %>% constructCredibleSet(thresh)}, 27, 0.99) %>%
  plyr::ldply(.id = NULL) %>% tbl_df()

