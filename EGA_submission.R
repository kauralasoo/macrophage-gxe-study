library("dplyr")
library("readr")
library("devtools")
load_all("macrophage-gxe-study/housekeeping/")

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







#Deal with acLDL samples
#Import complete line metadata
line_data = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") %>%
  dplyr::filter(donor != "mijn") %>%
  dplyr::select(-comment) %>%
  dplyr::select(-ips_culture_days, -passage_diff_bins) %>%
  dplyr::select(line_id, open_access) %>%
  unique()

acLDL_biosamples = readr::read_tsv("macrophage-gxe-study/data/sample_lists/biosamples/acLDL_sample_name_mapping.dates.with_accessions.txt", col_names = TRUE)
design_matrix = constructDesignMatrix_acLDL(acLDL_biosamples$sample_id) %>% tbl_df()
sample_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_acLDL_metadata.rds") %>%
  dplyr::select(line_id, donor)

is_open_access = dplyr::left_join(design_matrix, sample_metadata, by = "donor") %>% 
  dplyr::left_join(line_data, by = "line_id") %>%
  dplyr::select(sample_id, open_access)

acLDL_access_type = dplyr::left_join(acLDL_biosamples, is_open_access, by = "sample_id")

acLDL_open = dplyr::filter(acLDL_access_type, open_access == 1) %>% dplyr::select(-full_name)
acLDL_closed = dplyr::filter(acLDL_access_type, open_access == 0) %>% dplyr::select(-full_name)
write.table(acLDL_open, "results/sample_lists/submission/acLDL.open_access.txt", sep = ",", row.names = FALSE, quote = FALSE)
write.table(acLDL_closed, "results/sample_lists/submission/acLDL.managed_access.txt", sep = ",", row.names = FALSE, quote = FALSE)


#Import acLDL lanelets
acldl_lanelets = readr::read_tsv("macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_all.txt", col_names = c("sample_id", "lanelet_ids"))
acldl_lanelet_pairs = purrr::by_row(acldl_lanelets, 
                          ~data_frame(sample_id = .$sample_id, 
                          lanelet_id = unlist(strsplit(.$lanelet_ids, ";"))))$.out %>% 
  purrr::map_df(identity)

#Extract open and managed lanelets
open_lanelets = dplyr::semi_join(acldl_lanelet_pairs, acLDL_open, by = "sample_id") %>%
  dplyr::arrange(sample_id)
managed_lanelets = dplyr::semi_join(acldl_lanelet_pairs, acLDL_closed, by = "sample_id") %>%
  dplyr::arrange(sample_id)
write.table(open_lanelets, "results/sample_lists/submission/acLDL.open_access.lanelets.txt", sep = ",", row.names = FALSE, quote = FALSE)
write.table(managed_lanelets, "results/sample_lists/submission/acLDL.managed_access.lanelets.txt", sep = ",", row.names = FALSE, quote = FALSE)







