
conditions = c("naive","IFNg","SL1344","IFNg_SL1344")
chromosomes = as.character(c(1:22))

#Import all coloc results and merge them together
results = c()
for (condition in conditions){
  for (chr in chromosomes){
    coloc_file = paste("processed/ATAC_RNA/eQTL_caQTL_coloc/eQTL_caQTL_coloc", condition, chr, "txt", sep = ".")
    coloc_table = readr::read_tsv(coloc_file) %>%
      dplyr::mutate(condition_name = condition)
    results = dplyr::bind_rows(results, coloc_table)
  }
}

#Use reasonable filter to detect condfident colocs
filtered_hits = results %>%
  dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power) %>%
  dplyr::filter(PP_power > 0.8, PP_coloc > 0.9) %>%
  dplyr::group_by(gene_id, peak_id) %>%
  dplyr::arrange(gene_id, peak_id, -PP_coloc) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup()
saveRDS(filtered_hits, "results/ATAC_RNA_overlaps/QTL_overlap_list_coloc.rds")
