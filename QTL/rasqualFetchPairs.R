library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")

#Import SNP-peak pairs
rasqual_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalues_df = ldply(min_pvalue_hits, .id = "condition_name")
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Import SNP coords
snp_coords = readr::read_delim("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.snp_coords.txt", 
            delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))
selected_snp_coords = dplyr::filter(snp_coords, snp_id %in% unique(joint_pairs$snp_id)) %>% 
  dplyr::transmute(seqnames = chr, start = pos, end = pos, snp_id, strand = "*") %>%
  dataFrameToGRanges()
unique_genes = unique(joint_pairs$gene_id)

#Import SNP-peak pairs from disk
naive_snp_res = tabixFetchSNPs(selected_snp_coords, "results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz") %>% 
  dplyr::filter(gene_id %in% unique_genes) %>% 
  dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id"))
IFNg_snp_res = tabixFetchSNPs(selected_snp_coords, "results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz") %>% 
  dplyr::filter(gene_id %in% unique_genes) %>% 
  dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id"))
SL1344_snp_res = tabixFetchSNPs(selected_snp_coords, "results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz") %>% 
  dplyr::filter(gene_id %in% unique_genes) %>% 
  dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id"))
IFNg_SL1344_snp_res = tabixFetchSNPs(selected_snp_coords, "results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz") %>% 
  dplyr::filter(gene_id %in% unique_genes) %>% 
  dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id"))

rasqual_selected_results = list(naive = naive_snp_res, IFNg = IFNg_snp_res, SL1344 = SL1344_snp_res, IFNg_SL1344 = IFNg_SL1344_snp_res)
saveRDS(rasqual_selected_results, "results/ATAC/QTLs/rasqual_selected_results.rds")
