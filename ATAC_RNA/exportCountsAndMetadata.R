library("dplyr")

#Import complete line metadata
line_data = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") %>%
  dplyr::filter(donor != "mijn") %>%
  dplyr::select(-comment) %>%
  dplyr::select(-ips_culture_days, -passage_diff_bins)

#Import ATAC and RNA datasets
expression_list = readRDS("results/SL1344/combined_expression_data.rds")
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

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

#Export raw RNA-seq data counts
write.table(expression_list$counts, "figures/tables/RNA_count_matrix.txt", sep = "\t", 
            quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(expression_list$sample_metadata, "figures/tables/RNA_sample_metadata.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(expression_list$gene_metadata, "figures/tables/RNA_gene_metadata.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

#Export raw ATAC-seq data counts
write.table(atac_list$counts, "figures/tables/ATAC_count_matrix.txt", sep = "\t", 
            quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(atac_list$sample_metadata, "figures/tables/ATAC_sample_metadata.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(atac_list$gene_metadata, "figures/tables/ATAC_gene_metadata.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)


#Export ATAC-seq metadata before filtering
#Export sample metadata before filtering (to highlight QC criteria)
atac_biosamples = readr::read_tsv("macrophage-gxe-study/data/sample_lists/biosamples/atac_sample_name_mapping.dates.with_accessions.txt") %>%
  dplyr::transmute(sample_id, biosamples_accession = accession)
atac_qc_metadata = readRDS("macrophage-gxe-study/data/chromatin/ATAC/compiled_atac_metadata.rds") %>%
  dplyr::transmute(sample_id, donor, line_id, replicate, condition_name, atac_date = stimulation_date, 
                   assigned, submission_date, multiplex_pool, bioanalyzer, nextera_kit,
                   assigned_frac = round(assigned_frac,2), 
                   mt_frac = round(mt_frac,2), duplication = round(percent_duplication,2), 
                   peak_count, length_ratio = round(short_long_ratio,2)) %>%
  dplyr::mutate(qtl_mapping = ifelse(sample_id %in% c("fikt_B_ATAC", "qaqx_A_ATAC", "qaqx_B_ATAC", "uaqe_A_ATAC", "uaqe_B_ATAC"), "no","yes")) %>%
  dplyr::mutate(diff_access = ifelse(donor %in% 
                                       c("qolg","vass","kuxp","cicb", "febc","eiwy","oapg","nukw","hayt","bima","pamv","guss","eipl","iill","podx","pelm"),"yes","no")) %>%
  dplyr::arrange(sample_id) %>%
  dplyr::select(-donor) %>%
  dplyr::left_join(atac_biosamples, by = "sample_id")
write.table(atac_qc_metadata, "figures/tables/ATAC_QC_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

