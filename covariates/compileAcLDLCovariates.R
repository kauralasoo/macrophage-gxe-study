library("dplyr")
library("ggplot2")
library("tidyr")

#Import compiled line metadata
line_data = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")

#Import acLDL metadata
acLDL_dat = read.csv("macrophage-gxe-study/data/sample_lists/acLDL/acLDL_line_data.csv", stringsAsFactors = FALSE, na.strings = "") %>%
  dplyr::select(-comment, -library_pool, -macrophage_harvest) %>% tbl_df() %>%
  dplyr::mutate(acLDL_date = as.Date(acLDL_date, "%d/%m/%Y")) %>%
  dplyr::mutate(lalistat_date = as.Date(lalistat_date, "%d/%m/%Y")) %>%
  dplyr::mutate(rna_extraction = as.Date(rna_extraction, "%d/%m/%Y")) %>%
  dplyr::mutate(rna_submit = as.Date(rna_submit, "%d/%m/%Y"))

#Join tables together
line_data = dplyr::left_join(line_data, acLDL_dat, by = c("line_id", "replicate")) %>% 
  dplyr::filter(!is.na(acLDL_date)) %>%
  dplyr::arrange(line_id) %>%
  dplyr::select(-comment)

#Calculate duration of differentiation and make bins of 10 days
line_data = dplyr::mutate(line_data, macrophage_diff_days = as.numeric(acLDL_date - EB_formation)) %>%
  dplyr::mutate(macrophage_diff_bins = floor(macrophage_diff_days/10)*10) %>%
  dplyr::mutate(macrophage_diff_bins = ifelse(macrophage_diff_bins > 50, 50, macrophage_diff_bins)) #Cut-off at 50 days

#Remove flow purity for samples were flow was done more than 2 weeks after RNA extraction
line_data = dplyr::mutate(line_data, 
                          max_purity_filtered = ifelse(as.numeric(line_data$flow_date - line_data$acLDL_date) > 14, NA, max_purity),
                          mean_purity_filtered = ifelse(as.numeric(line_data$flow_date - line_data$acLDL_date) > 14,NA,mean_purity)) %>%
  dplyr::mutate(purity_bins = ifelse(max_purity_filtered < 0.98, "low", "high"))

#Save results to disk
saveRDS(line_data, "macrophage-gxe-study/data/covariates/compiled_acLDL_metadata.rds")
write.table(line_data, "macrophage-gxe-study/data/covariates/compiled_acLDL_metadata.txt", row.names = FALSE, sep = "\t", quote = FALSE)


