library("plyr")
library("dplyr")
library("purrr")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")

#Import SNP-peak pairs
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Import variant information
snp_coords = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
selected_snp_coords = dplyr::filter(snp_coords, snp_id %in% unique(joint_pairs$snp_id)) %>% 
  dplyr::transmute(seqnames = chr, start = pos, end = pos, snp_id, strand = "*") %>%
  dataFrameToGRanges()
unique_genes = unique(joint_pairs$gene_id)

#Import RASQUAL summary stats for SNP-peak pairs from disk
snp_res_list = purrr::map(qtlResults()$rna_rasqual, ~tabixFetchSNPs(selected_snp_coords, .) %>%
                            dplyr::filter(gene_id %in% unique_genes) %>%
                            dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id")))
saveRDS(snp_res_list, "results/SL1344/eQTLs/rasqual_selected_pvalues.rds")
