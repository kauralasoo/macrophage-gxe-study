library("devtools")
library("dplyr")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("../seqUtils/")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import eQTLs
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_qtl_df = extractQTLsFromList(rasqual_min_pvalues, fdr_cutoff = 0.1)
rasqual_pairs = dplyr::arrange(rasqual_qtl_df, p_fdr) %>%
  dplyr::select(gene_id, snp_id) %>% unique()
rasqual_filtered_pairs = filterHitsR2(rasqual_pairs, vcf_file$genotypes, .8)

#Import GWAS catalog
filtered_catalog = readRDS("../macrophage-gxe-study/annotations/gwas_catalog_v1.0.1-downloaded_2016-03-02.filtered.rds")

#GWAS overlaps only in the naive condition
naive_gwas_overlaps = dplyr::filter(rasqual_qtl_df, condition_name == "naive") %>% dplyr::select(gene_id, snp_id) %>% 
  unique() %>%
  findGWASOverlaps(filtered_catalog, vcf_file, min_r2 = 0.8)
naive_relative_olaps = rankTraitsByOverlapSize(naive_gwas_overlaps, filtered_catalog)
dplyr::filter(naive_relative_olaps, trait_size > 10, overlap_size >= 5) %>% arrange(-fraction)

#Find overlaps between lead variants and eQTLs
all_gwas_overlaps = dplyr::select(rasqual_filtered_pairs, gene_id, snp_id) %>% unique() %>%
  findGWASOverlaps(filtered_catalog, vcf_file, min_r2 = 0.8)

relative_olaps = rankTraitsByOverlapSize(all_gwas_overlaps, filtered_catalog)
dplyr::filter(relative_olaps, trait_size > 10, overlap_size >= 5) %>% arrange(-fraction)

#Find overlaps with splicing QTLs
#Import Salmon QTLs
salmon_qtl_hits = readRDS("results/SL1344/salmon/salmon_qtl_hits.rds")
salmon_qtl_df = purrr::map_df(salmon_qtl_hits, identity, .id = "condition_name")
salmon_pairs = dplyr::arrange(salmon_qtl_df, p_fdr) %>%
  dplyr::transmute(gene_id = ensembl_gene_id, snp_id) %>% 
  unique()
salmon_filtered_pairs = filterHitsR2(salmon_pairs, vcf_file$genotypes, .8)

#Salmon naive overlaps
salmon_naive_overlaps = dplyr::filter(salmon_qtl_df, condition_name == "naive") %>% 
  dplyr::transmute(gene_id = ensembl_gene_id, snp_id) %>% 
  unique() %>%
  findGWASOverlaps(filtered_catalog, vcf_file, min_r2 = 0.8)
salmon_naive_olaps = rankTraitsByOverlapSize(salmon_naive_overlaps, filtered_catalog)
dplyr::filter(salmon_naive_olaps, trait_size > 10, overlap_size >= 5) %>% arrange(-fraction)x

#Find overlaps between lead variants and eQTLs
salmon_gwas_overlaps = findGWASOverlaps(salmon_filtered_pairs, filtered_catalog, vcf_file, min_r2 = 0.8)
relative_olaps = rankTraitsByOverlapSize(salmon_gwas_overlaps, filtered_catalog)
dplyr::filter(relative_olaps, trait_size > 10, overlap_size >= 5) %>% arrange(-fraction)

#All expression QTLs
all_pairs = dplyr::bind_rows(dplyr::transmute(salmon_qtl_df, gene_id = ensembl_gene_id, snp_id, p_fdr) , 
                 dplyr::select(rasqual_qtl_df, gene_id, snp_id, p_fdr)) %>% 
  group_by(gene_id) %>% 
  arrange(gene_id, p_fdr) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(gene_id, snp_id) %>% 
  unique()
all_filtered_pairs = filterHitsR2(all_pairs, vcf_file$genotypes, .8)

#Rank traits by relative overlap
combined_gwas_overlaps = findGWASOverlaps(all_filtered_pairs, filtered_catalog, vcf_file, min_r2 = 0.8)
relative_olaps = rankTraitsByOverlapSize(combined_gwas_overlaps, filtered_catalog)
ranked_overlaps = dplyr::filter(relative_olaps, trait_size > 15) %>% 
  arrange(-fraction) %>% 
  dplyr::mutate(fraction = round(fraction,2))
write.table(ranked_overlaps, "figures/supplementary/GWAS_ranked_traits.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#Count GWAS overlaps from different analyses
overlap_counts = data_frame(condition = c("RASQUAL (naive)","RASQUAL (all)", "Salmon (naive)", 
                                          "Salmon (all)", "RASQUAL + Salmon"),
           overlap_count = c(nrow(naive_gwas_overlaps), nrow(all_gwas_overlaps), nrow(salmon_naive_olaps), 
                             nrow(salmon_gwas_overlaps), nrow(combined_gwas_overlaps))) %>%
  dplyr::mutate(condition = factor(condition, levels = condition))

overlap_count_plot = ggplot(overlap_counts, aes(x = condition, y = overlap_count)) + 
  geom_bar(stat = "identity") + 
  theme_light() + 
  ylab("Number of overlaps") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 15))
ggsave("figures/supplementary/GWAS_overap_counts.pdf", plot = overlap_count_plot, width = 6, height = 4.5)


#IRF5 overlaps
irf5_variants = data_frame(gene_id = "IRF5", snp_id = c("rs199508964","rs3778754","rs10954213"))
irf5_overlaps = findGWASOverlaps(irf5_variants, filtered_catalog, vcf_file, min_r2 = 0.8)

