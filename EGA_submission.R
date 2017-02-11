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


#Export all missing lanelets for open access and managed access samples
rna_lanelets = readr::read_tsv("macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt", col_names = c("sample_id", "lanelet_ids"))
rna_lanelet_pairs = purrr::by_row(rna_lanelets, 
    ~data_frame(sample_id = .$sample_id, 
   lanelet_id = unlist(strsplit(.$lanelet_ids, ";"))))$.out %>% 
  purrr::map_df(identity)

atac_lanelets = readr::read_tsv("macrophage-gxe-study/data/chromatin/ATAC/ATAC_Salmonella_names.txt", col_names = c("sample_id", "lanelet_ids"))
atac_lanelet_pairs = purrr::by_row(atac_lanelets, 
      ~data_frame(sample_id = .$sample_id, 
      lanelet_id = unlist(strsplit(.$lanelet_ids, ";"))))$.out %>% 
  purrr::map_df(identity)


#Import open access samples
rna_open = readr::read_tsv("results/sample_lists/submission/salmonella_rna.open_access.txt") %>%
  dplyr::select(-lanelet_id) %>%
  dplyr::left_join(rna_lanelet_pairs, by = "sample_id")
atac_open = readr::read_tsv("results/sample_lists/submission/salmonella_atac.open_access.txt") %>%
  dplyr::select(-lanelet_id) %>%
  dplyr::left_join(atac_lanelet_pairs, by = "sample_id")

write.table(rna_open, "results/sample_lists/submission/salmonella_rna.open_access.all_lanelets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(atac_open, "results/sample_lists/submission/salmonella_atav.open_access.all_lanelets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#Import unknown lanelets
unknown_lanelets = readr::read_tsv("results/sample_lists/4584_files.txt", col_names = "cram_file") %>%
  tidyr::separate(cram_file, c("lanelet_id", "suffix"), sep = "\\.")

#Export sample names and accessions for missing lanelets
rna_matched = dplyr::semi_join(rna_open, unknown_lanelets, by = "lanelet_id") %>% dplyr::arrange(sample_id)
atac_matched = dplyr::semi_join(atac_open, unknown_lanelets, by = "lanelet_id") %>% dplyr::arrange(sample_id)

write.table(rna_matched, "results/sample_lists/submission/salmonella_rna.open_access.missing_lanelets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(atac_matched, "results/sample_lists/submission/salmonella_atav.open_access.missing_lanelets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


#Import managed access samples
rna_managed = readr::read_tsv("results/sample_lists/submission/salmonella_rna.managed_access.txt") %>%
  dplyr::select(-lanelet_id) %>%
  dplyr::left_join(rna_lanelet_pairs, by = "sample_id")
atac_managed = readr::read_tsv("results/sample_lists/submission/salmonella_atac.managed_access.txt") %>%
  dplyr::select(-lanelet_id) %>%
  dplyr::left_join(atac_lanelet_pairs, by = "sample_id")


write.table(rna_managed, "results/sample_lists/submission/salmonella_rna.managed_access.all_lanelets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(atac_managed, "results/sample_lists/submission/salmonella_atac.managed_access.all_lanelets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


