library("dplyr")
library("tidyr")

#Import complete line metadata
line_data = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") %>%
  dplyr::filter(donor != "mijn") %>%
  dplyr::select(open_access, donor) %>% unique()

#Import biosamples lists
sl1344_biosamples = readr::read_tsv("macrophage-gxe-study/data/sample_lists/biosamples/salmonella_sample_name_mapping.dates.with_accessions.txt") %>%
  dplyr::select(-full_name)
atac_biosamples = readr::read_tsv("macrophage-gxe-study/data/sample_lists/biosamples/atac_sample_name_mapping.dates.with_accessions.txt") %>%
  dplyr::select(-full_name)

#Add open access info for RNA
sl1344_new = tidyr::separate(sl1344_biosamples, sample_id, c("donor", "clone"), sep = "_", remove = FALSE) %>% 
  dplyr::left_join(line_data, by = "donor") %>% dplyr::select(-donor, -clone)
write.table(dplyr::filter(sl1344_new, open_access == 1), 
            "results/sample_lists/submission/salmonella_rna.open_access.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dplyr::filter(sl1344_new, open_access == 0), 
            "results/sample_lists/submission/salmonella_rna.managed_access.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

#Add open access info for ATAC
atac_new = tidyr::separate(atac_biosamples, sample_id, c("donor", "clone"), sep = "_", remove = FALSE) %>% 
  dplyr::left_join(line_data, by = "donor") %>% dplyr::select(-donor, -clone)
write.table(dplyr::filter(atac_new, open_access == 1), 
            "results/sample_lists/submission/salmonella_atac.open_access.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dplyr::filter(atac_new, open_access == 0), 
            "results/sample_lists/submission/salmonella_atac.managed_access.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

