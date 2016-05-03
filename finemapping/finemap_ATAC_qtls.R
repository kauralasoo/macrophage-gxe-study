library("devtools")
library("purrr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")


#List of ATAC tabix files
atac_tabix_list = list(naive = "../macrophage-chromatin/results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz",
                       IFNg = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz",
                       SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz",
                       IFNg_SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz")

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Fetch credible sets for each peak in each condition
naive_cs = fetchCredibleSets(min_pvalues_hits$naive[1:10,], atac_list$gene_metadata, atac_list$gene_metadata, atac_tabix_list$naive, vcf_file, cis_window = 5e4)
IFNg_cs = fetchCredibleSets(min_pvalues_hits$IFNg[1:10,], atac_list$gene_metadata, atac_list$gene_metadata, atac_tabix_list$IFNg, vcf_file, cis_window = 5e4)
SL1344_cs = fetchCredibleSets(min_pvalues_hits$SL1344[1:10,], atac_list$gene_metadata, atac_list$gene_metadata, atac_tabix_list$SL1344, vcf_file, cis_window = 5e4)
IFNg_SL1344_cs = fetchCredibleSets(min_pvalues_hits$IFNg_SL1344[1:10,], atac_list$gene_metadata, atac_list$gene_metadata, atac_tabix_list$IFNg_SL1344, vcf_file, cis_window = 5e4)

#Save credible sets to disk
credible_set_list = list(naive = naive_cs, IFNg = IFNg_cs, SL1344 = SL1344_cs, IFNg_SL1344 = IFNg_SL1344_cs)
saveRDS(credible_set_list, "results/ATAC/QTLs/rasqual_credible_sets.rds")

