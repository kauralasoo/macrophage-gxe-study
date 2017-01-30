library("dplyr")

#Import complete line metadata
line_data = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") %>%
  dplyr::filter(donor != "mijn") %>%
  dplyr::select(-comment) %>%
  dplyr::select(-ips_culture_days, -passage_diff_bins) %>%
  dplyr::select(donor, genotype_id, sex) %>% unique()

#Process RNA data
ega_meta = readr::read_tsv("results/sample_lists/submission/salmonella_rna.managed_access.txt")
rna_design = constructDesignMatrix_SL1344(ega_meta$sample_id) %>% dplyr::select(sample_id, donor) %>% tbl_df()

#Construct additional metadata
rna_additional_meta = dplyr::left_join(ega_meta, rna_design, by = "sample_id") %>% 
  dplyr::left_join(line_data, by = "donor") %>% 
  dplyr::transmute(accession, donor = genotype_id, phenotype = "RNA-seq", gender = sex)
write.table(rna_additional_meta, "results/sample_lists/submission/salmonella_rna.additional_metdata.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#Process ATAC data
atac_meta = readr::read_tsv("results/sample_lists/submission/salmonella_atac.managed_access.txt")
atac_design = constructDesignMatrix_ATAC(atac_meta$sample_id) %>% dplyr::select(sample_id, donor) %>% tbl_df()

#Construct additional metadata
atac_additional_meta = dplyr::left_join(atac_meta, atac_design, by = "sample_id") %>% 
  dplyr::left_join(line_data, by = "donor") %>% 
  dplyr::transmute(accession, donor = genotype_id, phenotype = "ATAC-seq", gender = sex)
write.table(atac_additional_meta, "results/sample_lists/submission/salmonella_atac.additional_metdata.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#Import shared samples
shared_samples = readr::read_tsv("results/sample_lists/shared_samples.txt", col_names = "accession")
rna_lanelets = readr::read_tsv("macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt", col_names = c("sample_id", "lanelet_ids"))
atac_lanelets = readr::read_tsv("macrophage-gxe-study/data/chromatin/ATAC/ATAC_Salmonella_names.txt", col_names = c("sample_id", "lanelet_ids"))

rna_shared_lanelets = readr::read_tsv("macrophage-gxe-study/data/sample_lists/biosamples/salmonella_sample_name_mapping.dates.with_accessions.txt") %>% 
  dplyr::select(sample_id, accession) %>% 
  dplyr::semi_join(shared_samples, by = "accession") %>%
  dplyr::left_join(rna_lanelets, by = "sample_id")
atac_shared_lanelets = readr::read_tsv("macrophage-gxe-study/data/sample_lists/biosamples/atac_sample_name_mapping.dates.with_accessions.txt") %>% 
  dplyr::select(sample_id, accession) %>% 
  dplyr::semi_join(shared_samples, by = "accession") %>%
  dplyr::left_join(atac_lanelets, by = "sample_id")

#Import full lanelet lists
write.table(rna_shared_lanelets, "results/sample_lists/submission/rna_shared_lanelets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(atac_shared_lanelets, "results/sample_lists/submission/atac_shared_lanelets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


