library("dplyr")

#Import complete line metadata
line_data = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") %>%
  dplyr::filter(donor != "mijn") %>%
  dplyr::select(-comment) %>%
  dplyr::select(-ips_culture_days, -passage_diff_bins)

#Import ATAC and RNA datasets
expression_list = readRDS("results/SL1344/combined_expression_data_covariates.rds")
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Extract lines with good RNA-seq data
rna_lines = dplyr::select(expression_list$sample_metadata, line_id, replicate) %>% 
  unique() %>% 
  dplyr::mutate(RNA_included = 1)

#Extract lines with good ATAC-seq data
atac_lines = dplyr::select(atac_list$sample_metadata, line_id, replicate) %>% 
  unique() %>% 
  dplyr::mutate(ATAC_included = 1)

#Add RNA and ATAC flags to the metadata:
compiled_data = dplyr::left_join(line_data, rna_lines, by = c("line_id","replicate")) %>% 
  dplyr::left_join(atac_lines, by = c("line_id","replicate")) %>% 
  dplyr::mutate(RNA_included = ifelse(is.na(RNA_included), 0, 1)) %>%
  dplyr::mutate(ATAC_included = ifelse(is.na(ATAC_included), 0, 1))

#Export line metadata
write.table(compiled_data, "figures/tables/line_differentiation_metadata.txt", 
           sep = "\t", quote = FALSE, row.names = FALSE)