
#ADD RNA conentrations
rna_concentrations = readRDS("macrophage-gxe-study/data/covariates/rna_concentrations.rds")


#Add mean RNA concentration and divide it into bins (of 100 ng/ul)
mean_rna_concentrations = dplyr::select(rna_concentrations, donor, ng_ul_mean, replicate) %>% unique()
rna_bins = data_frame(rna_bin = c(0,1,2), rna_concentration = c("0-100","100-200","200+"))
rna_cons = dplyr::mutate(mean_rna_concentrations, rna_bin = floor(mean_rna_concentrations$ng_ul_mean/100)) %>% 
  dplyr::mutate(rna_bin = ifelse(rna_bin > 2, 2, rna_bin)) %>% #Remove anything bigger than 300
  dplyr::left_join(rna_bins, by = "rna_bin") %>% 
  dplyr::select(donor, replicate, rna_concentration, ng_ul_mean)





#Add iPS cell culture duration
line_data = 
  dplyr::mutate(mf_diff_days = as.numeric(salmonella - MF_harvest)) #MF diff duration




#Calculate duration of differentiation and make bins of 10 days
line_data = dplyr::mutate(line_data, diff_days = as.numeric(salmonella - EB_formation)) %>%
  dplyr::mutate(diff_bins = floor(diff_days/10)*10) %>%
  dplyr::mutate(diff_bins = ifelse(diff_bins > 50, 50, diff_bins)) #Cut-off at 50 days



#Remove flow purity for samples were flow was done more than 2 weeks after RNA extraction
line_data = dplyr::mutate(line_data, max_purity_filtered = ifelse(as.numeric(line_data$flow_date - line_data$salmonella) > 14, NA, max_purity),
                          mean_purity_filtered = ifelse(as.numeric(line_data$flow_date - line_data$salmonella) > 14,NA,mean_purity)) %>%
  dplyr::mutate(purity_bins = ifelse(max_purity < 0.98, "low", "high"))


