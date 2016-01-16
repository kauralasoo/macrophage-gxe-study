line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")
genotype_id_list = dplyr::filter(line_metadata, !is.na(salmonella_date) | !is.na(atac_date)) %>% 
  dplyr::select(line_id, genotype_id, salmonella_date, atac_date) %>% 
  dplyr::filter(line_id != "fpdj_3", line_id != "mijn_2") %>% 
  dplyr::arrange(line_id) %>%
  dplyr::select(genotype_id) %>% unique()

write.table(genotype_id_list$genotype_id, "macrophage-gxe-study/data/sample_lists/SL1344/SL1344_gt_list.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)