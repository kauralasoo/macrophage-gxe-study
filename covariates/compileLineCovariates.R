library("dplyr")
library("ggplot2")
library("tidyr")

#Import different files
dat = tbl_df(read.csv("macrophage-gxe-study/data/sample_lists/line_metadata.csv", stringsAsFactors = FALSE, na.strings = ""))
flow_purity = readRDS("macrophage-gxe-study/data/covariates/flow_cytometry_purity.rds")
gender_map = read.table("macrophage-gxe-study/data/sample_lists/line_gender_map.txt", header = TRUE, stringsAsFactors = FALSE)
genotypes = read.table("macrophage-gxe-study/data/sample_lists/genotype_sample_names.txt", stringsAsFactors = FALSE)
isOpenAccess = read.table("macrophage-gxe-study/data/sample_lists/HipSci_is_open_access.txt", stringsAsFactors = FALSE)

#Convert all date fields to date
line_metadata = dat %>%
  tidyr::separate(line_id, into = c("donor","clone"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(ips_received = as.Date(ips_received, "%d/%m/%Y"), 
                ips_started = as.Date(ips_started, "%d/%m/%Y"),
                EB_formation = as.Date(EB_formation, "%d/%m/%Y"),
                diff_start = as.Date(diff_start, "%d/%m/%Y"),
                salmonella_date = as.Date(salmonella_date, "%d/%m/%Y"),
                atac_date = as.Date(atac_date, "%d/%m/%Y"),
                flow_date = as.Date(flow_date, "%d/%m/%Y"),
                terminated = as.Date(terminated, "%d/%m/%Y"))

#Add flow purity
mean_flow_purity = flow_purity %>%
  dplyr::filter(!(donor == "gedo" & channel == "Pacific.Blue.A")) %>% #Remove an outlier measurement
  group_by(donor,flow_date) %>% 
  dplyr::summarize(mean_purity = mean(purity), max_purity = max(purity))

#Add gender to the metadata
tidy_sex = dplyr::mutate(gender_map, new_gender = ifelse(gender == "Female" | gender == "female*", "female", gender)) %>% 
  dplyr::mutate(new_gender = ifelse(gender == "Male" | gender == "male*", "male", new_gender)) %>% 
  dplyr::select(line_id, new_gender) %>% 
  dplyr::rename(sex = new_gender) %>%
  unique()

#Add genotype ids and open access status
colnames(isOpenAccess) = c("genotype_id", "open_access")
genotypes_open_access_names = tidyr::separate(isOpenAccess, genotype_id, into = c("batch_id", "line_id"), sep = "-", remove = FALSE) %>% 
  dplyr::select(genotype_id, line_id, open_access)

#Put all of the annotations together
line_data = dplyr::left_join(line_metadata, mean_flow_purity, by = c("donor", "flow_date")) %>%
  dplyr::left_join(genotypes_open_access_names, by = "line_id") %>%
  dplyr::left_join(tidy_sex, by = "line_id") %>%
  dplyr::mutate(ips_culture_days = as.numeric(EB_formation - ips_started)) #iPS culture duration

#Split the iPS passage into 3 bins
line_data = dplyr::mutate(line_data, passage_diff_bins = floor(line_data$passage_diff/10)*10) %>% 
  dplyr::mutate(passage_diff_bins = ifelse(passage_diff_bins < 20, 20, passage_diff_bins)) %>% #Below 30 is lowest
  dplyr::mutate(passage_diff_bins = ifelse(passage_diff_bins > 40, 40, passage_diff_bins)) #Above 40 is highest

#Save line metadata to disk
saveRDS(line_data, "macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")
write.table(line_data, "macrophage-gxe-study/data/covariates/compiled_line_metadata.txt", row.names = FALSE, sep = "\t", quote = FALSE)

