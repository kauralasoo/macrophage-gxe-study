library("dplyr")
library("readr")
library("devtools")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Load all data

#Import sample lists
sl1344_samples = readr::read_delim("macrophage-gxe-study/data/sample_lists/SL1344/sample_name_mapping.release.txt", delim = "\t", col_names = TRUE)
acLDL_samples = readr::read_delim("macrophage-gxe-study/data/sample_lists/acLDL/acLDL_sample_name_mapping.release.txt", delim = "\t", col_names = TRUE)
atac_samples = readr::read_delim("macrophage-gxe-study/data/chromatin/ATAC/ATAC_sample_name_mapping.txt", delim = "\t", col_names = TRUE)

#Construct design matrixt
sl1344_design = constructDesignMatrix_SL1344(sl1344_samples$sample_id) %>%
  dplyr::mutate(replicate = ifelse(donor == "babk",2,replicate))#Change babk replicate to two
acLDL_design = constructDesignMatrix_acLDL(acLDL_samples$sample_id)
atac_design = constructDesignMatrix_ATAC(atac_samples$sample_id)

#Import metadata
sl1344_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_salmonella_metadata.rds")
acLDL_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_acLDL_metadata.rds")
atac_metadata = readRDS("macrophage-gxe-study/data/chromatin/ATAC/compiled_atac_metadata.rds")

#Add dates
sl1344_dates = dplyr::left_join(sl1344_samples, sl1344_design, by = "sample_id") %>% 
  dplyr::left_join(sl1344_metadata, by = c("donor","replicate")) %>% 
  dplyr::select(sample_id, lanelet_id, sequencescape_id, sanger_id, diff_start, salmonella_date)
acLDL_dates = dplyr::left_join(acLDL_samples, acLDL_design, by = "sample_id") %>% 
  dplyr::left_join(acLDL_metadata, by = c("donor")) %>% 
  dplyr::select(sample_id, lanelet_id, sequencescape_id, sanger_id, diff_start, acLDL_date)
atac_dates = dplyr::left_join(atac_samples, atac_design, by = "sample_id") %>% 
  dplyr::left_join(atac_metadata, by = c("sample_id")) %>% 
  dplyr::select(sample_id, lanelet_id, sequencescape_id, sanger_id, diff_start, atac_date)

#Save dates
write.table(acLDL_dates, "results/sample_lists/acLDL_sample_name_mapping.dates.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(atac_dates, "results/sample_lists/atac_sample_name_mapping.dates.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sl1344_dates, "results/sample_lists/salmonella_sample_name_mapping.dates.txt", sep = "\t", quote = FALSE, row.names = FALSE)


