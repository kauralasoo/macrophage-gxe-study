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
min_pvalue_hits = lapply(min_pvalue_list, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = ldply(min_pvalue_hits, .id = "condition_name") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(-chisq)
joint_pairs = dplyr::select(min_pvalue_df, gene_id, snp_id) %>% unique()
saveRDS(min_pvalue_list, "results/SL1344/eQTLs/rasqual_min_pvalues.rds")

#Import SNP coords
snp_coords = readr::read_delim("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz", 
                               delim = "\t", col_types = "cdccc", col_names = c("chr","pos","snp_id","ref","alt"))
selected_snp_coords = dplyr::filter(snp_coords, snp_id %in% unique(joint_pairs$snp_id)) %>% 
  dplyr::transmute(seqnames = chr, start = pos, end = pos, snp_id, strand = "*") %>%
  dataFrameToGRanges()
unique_genes = unique(joint_pairs$gene_id)

#Import SNP-peak pairs from disk
naive_snp_res = tabixFetchSNPs(selected_snp_coords, "results/SL1344/rasqual/output/naive_500kb/naive_500kb.sorted.txt.gz") %>% 
  dplyr::filter(gene_id %in% unique_genes) %>% 
  dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id"))
IFNg_snp_res = tabixFetchSNPs(selected_snp_coords, "results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.sorted.txt.gz") %>% 
  dplyr::filter(gene_id %in% unique_genes) %>% 
  dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id"))
SL1344_snp_res = tabixFetchSNPs(selected_snp_coords, "results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.sorted.txt.gz") %>% 
  dplyr::filter(gene_id %in% unique_genes) %>% 
  dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id"))
IFNg_SL1344_snp_res = tabixFetchSNPs(selected_snp_coords, "results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.sorted.txt.gz") %>% 
  dplyr::filter(gene_id %in% unique_genes) %>% 
  dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id"))

rasqual_selected_results = list(naive = naive_snp_res, IFNg = IFNg_snp_res, SL1344 = SL1344_snp_res, IFNg_SL1344 = IFNg_SL1344_snp_res)
saveRDS(rasqual_selected_results, "results/SL1344/eQTLs/rasqual_selected_pvalues.rds")
