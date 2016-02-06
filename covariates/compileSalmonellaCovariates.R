library("dplyr")
library("ggplot2")
library("tidyr")

#Import compiled line metadata
line_data = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")

#Import Salmonella metadata
salmonella_dat = read.csv("macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_metadata.csv", 
                          stringsAsFactors = FALSE, na.strings = "") %>%
  tbl_df() %>%
  dplyr::mutate(macrophage_harvest = as.Date(macrophage_harvest, "%d/%m/%Y")) %>%
  dplyr::mutate(salmonella_date = as.Date(salmonella_date, "%d/%m/%Y")) %>%
  dplyr::mutate(rna_extraction = as.Date(rna_extraction, "%d/%m/%Y")) %>%
  dplyr::mutate(rna_submit = as.Date(rna_submit, "%d/%m/%Y")) %>%
  dplyr::mutate(rna_auto = ifelse(chemistry == "V4_auto", TRUE, FALSE))

#ADD RNA conentrations
rna_concentrations = readRDS("macrophage-gxe-study/data/covariates/rna_concentrations.rds")

#Add mean RNA concentration and divide it into bins (of 100 ng/ul)
mean_rna_concentrations = dplyr::select(rna_concentrations, donor, ng_ul_mean, replicate) %>% unique()
rna_bins = data_frame(rna_bin = c(0,1,2), rna_concentration = c("0-100","100-200","200+"))
rna_cons = dplyr::mutate(mean_rna_concentrations, rna_bin = floor(mean_rna_concentrations$ng_ul_mean/100)) %>% 
  dplyr::mutate(rna_bin = ifelse(rna_bin > 2, 2, rna_bin)) %>% #Remove anything bigger than 300
  dplyr::left_join(rna_bins, by = "rna_bin") %>% 
  dplyr::select(donor, replicate, rna_concentration, ng_ul_mean)

#Add to the line data
line_data = dplyr::left_join(line_data, rna_cons, by = c("donor","replicate")) %>%
  dplyr::left_join(salmonella_dat, by = c("line_id", "replicate", "salmonella_date")) %>%
  dplyr::select(-comment.x, -comment.y) %>%
  dplyr::filter(!is.na(salmonella_date)) %>%
  dplyr::filter(status == "Success")

#Calculate duration of differentiation and make bins of 10 days
line_data = dplyr::mutate(line_data, macrophage_diff_days = as.numeric(salmonella_date - EB_formation)) %>%
  dplyr::mutate(macrophage_diff_bins = floor(macrophage_diff_days/10)*10) %>%
  dplyr::mutate(macrophage_diff_bins = ifelse(macrophage_diff_bins > 50, 50, macrophage_diff_bins)) #Cut-off at 50 days

#Add iPS cell culture duration
line_data = dplyr::mutate(line_data, harvest_stimulation_days = as.numeric(salmonella_date - macrophage_harvest)) #MF diff duration

#Remove flow purity for samples were flow was done more than 2 weeks after RNA extraction
line_data = dplyr::mutate(line_data, 
              max_purity_filtered = ifelse(as.numeric(line_data$flow_date - line_data$salmonella_date) > 14, NA, max_purity),
              mean_purity_filtered = ifelse(as.numeric(line_data$flow_date - line_data$salmonella_date) > 14,NA,mean_purity)) %>%
  dplyr::mutate(purity_bins = ifelse(max_purity < 0.98, "low", "high"))

#Save results to disk
saveRDS(line_data, "macrophage-gxe-study/data/covariates/compiled_salmonella_metadata.rds")
write.table(line_data, "macrophage-gxe-study/data/covariates/compiled_salmonella_metadata.txt", row.names = FALSE, sep = "\t", quote = FALSE)

