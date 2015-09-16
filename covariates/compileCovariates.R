library("dplyr")
library("ggplot2")
library("tidyr")

#Import different files
dat = tbl_df(read.csv("macrophage-gxe-study/data/sample_lists/line_metadata.csv", stringsAsFactors = FALSE, na.strings = ""))
flow_purity = readRDS("macrophage-gxe-study/data/covariates/flow_cytometry_purity.rds")
rna_concentrations = readRDS("macrophage-gxe-study/data/covariates/rna_concentrations.rds")
gender_map = read.table("macrophage-gxe-study/data/sample_lists/line_gender_map.txt", header = TRUE, stringsAsFactors = FALSE)
genotypes = read.table("macrophage-gxe-study/data/sample_lists/genotype_sample_names.txt", stringsAsFactors = FALSE)

#Load genotype sample names
colnames(genotypes) = c("genotype_id")
genotypes_names = tidyr::separate(genotypes, genotype_id, into = c("batch_id", "line_id"), sep = "-", remove = FALSE) %>% 
  dplyr::select(genotype_id, line_id)

#Convert all date fields to date
line_metadata = dat %>%
  tidyr::separate(line_id, into = c("donor","clone"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(ips_received = as.Date(ips_received, "%d/%m/%Y"), 
                ips_started = as.Date(ips_started, "%d/%m/%Y"),
                EB_formation = as.Date(EB_formation, "%d/%m/%Y"),
                diff_start = as.Date(diff_start, "%d/%m/%Y"),
                MF_harvest = as.Date(MF_harvest, "%d/%m/%Y"),
                salmonella = as.Date(salmonella, "%d/%m/%Y"),
                flow_date = as.Date(flow_date, "%d/%m/%Y"),
                terminated = as.Date(terminated, "%d/%m/%Y"),
                extraction_date = as.Date(extraction_date, "%d/%m/%Y"),
                rna_submit = as.Date(rna_submit, "%d/%m/%Y"))

#Add mean RNA concentration and divide it into bins (of 100 ng/ul)
mean_rna_concentrations = dplyr::select(rna_concentrations, donor, ng_ul_mean, replicate) %>% unique()
rna_bins = data_frame(rna_bin = c(0,1,2), rna_concentration = c("0-100","100-200","200+"))
rna_cons = dplyr::mutate(mean_rna_concentrations, rna_bin = floor(mean_rna_concentrations$ng_ul_mean/100)) %>% 
  dplyr::mutate(rna_bin = ifelse(rna_bin > 2, 2, rna_bin)) %>% #Remove anything bigger than 300
  dplyr::left_join(rna_bins, by = "rna_bin") %>% 
  dplyr::select(donor, replicate, rna_concentration, ng_ul_mean)

#Add flow purity
mean_flow_purity = flow_purity %>%
  dplyr::filter(!(donor == "gedo" & channel == "Pacific.Blue.A")) %>% #Remove an outlier measurement
  group_by(donor,flow_date) %>% 
  dplyr::summarize(mean_purity = mean(purity), max_purity = max(purity))

#Add gender to the metadata
tidy_gender = dplyr::mutate(gender_map, new_gender = ifelse(gender == "Female" | gender == "female*", "female", gender)) %>% 
  dplyr::mutate(new_gender = ifelse(gender == "Male" | gender == "male*", "male", new_gender)) %>% 
  dplyr::select(line_id, new_gender) %>% 
  dplyr::rename(gender = new_gender) %>%
  unique()

#Put all of the annotations together
line_data = dplyr::left_join(line_metadata, rna_cons, by = c("donor", "replicate")) %>%
  dplyr::left_join(mean_flow_purity, by = c("donor", "flow_date")) %>%
  dplyr::left_join(genotypes_names, by = "line_id") %>%
  dplyr::left_join(tidy_gender, by = "line_id")

#Calculate duration of differentiation and make bins of 10 days
line_data = dplyr::mutate(line_data, diff_days = as.numeric(salmonella - EB_formation)) %>%
  dplyr::mutate(diff_bins = floor(diff_days/10)*10) %>%
  dplyr::mutate(diff_bins = ifelse(diff_bins > 50, 50, diff_bins)) #Cut-off at 50 days

#Add iPS cell culture duration
line_data = dplyr::mutate(line_data, ips_culture_days = as.numeric(EB_formation - ips_started)) %>% #iPS culture duration
  dplyr::mutate(mf_diff_days = as.numeric(salmonella - MF_harvest)) #MF diff duration

#Split the iPS passage into 3 bins
line_data = dplyr::mutate(line_data, passage_diff_bins = floor(line_data$passage_diff/10)*10) %>% 
  dplyr::mutate(passage_diff_bins = ifelse(passage_diff_bins < 20, 20, passage_diff_bins)) %>% #Below 30 is lowest
  dplyr::mutate(passage_diff_bins = ifelse(passage_diff_bins > 40, 40, passage_diff_bins)) #Above 40 is highest

#Remove flow purity for samples were flow was done more than 2 weeks after RNA extraction
line_data = dplyr::mutate(line_data, max_purity_filtered = ifelse(as.numeric(line_data$flow_date - line_data$salmonella) > 14, NA, max_purity),
                          mean_purity_filtered = ifelse(as.numeric(line_data$flow_date - line_data$salmonella) > 14,NA,mean_purity)) %>%
  dplyr::mutate(purity_bins = ifelse(max_purity < 0.98, "low", "high"))

saveRDS(line_data, "macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")
write.table(line_data, "macrophage-gxe-study/data/covariates/compiled_line_metadata.txt", row.names = FALSE, sep = "\t", quote = FALSE)

