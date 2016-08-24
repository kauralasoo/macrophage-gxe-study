#Import line metadata from disk
line_data = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") %>%
  dplyr::filter(!(donor %in% c("mijn")))

#Extract acLDL samples
acldl_success_samples = readRDS("macrophage-gxe-study/data/covariates/compiled_acLDL_metadata.rds") %>%
  dplyr::filter(!(donor %in% c("mijn", "xegx"))) %>%
  dplyr::select(line_id, replicate, acLDL_date)

#Add acldl data to line data
full_data = dplyr::left_join(line_data, acldl_success_samples, by = c("line_id","replicate")) %>%
  #If acLDL data present then change status to success
  dplyr::mutate(status = ifelse(!is.na(acLDL_date), "Success", status)) %>%
  dplyr::mutate(diff_fail = ifelse(status == "Diff_fail", 1, 0)) %>%
  dplyr::mutate(assay_success = ifelse(status == "Success", 1, 0)) %>%
  dplyr::mutate(flow_success = ifelse(is.na(flow_date), 0, 1))

#Count the number of differentiation attempts and the number of failures
diff_fail_count = dplyr::group_by(full_data, line_id) %>% 
  dplyr::summarise(diff_attempt_count = length(line_id), 
            diff_fail_count = sum(diff_fail), 
            assay_success_count = sum(assay_success),
            flow_success_count = sum(flow_success)) %>% 
  dplyr::arrange(-diff_attempt_count)
write.table(diff_fail_count, "macrophage-gxe-study/data/covariates/macrophage_diff_fail_count.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE)

#Make some queries
dplyr::filter(diff_fail_count, diff_attempt_count > diff_fail_count)
dplyr::filter(diff_fail_count, assay_success_count > 0)
dplyr::filter(diff_fail_count, flow_success_count > 0)

#Count the number of feeder and feeder free lines
media_types = dplyr::select(line_data, line_id, media) %>% 
  dplyr::group_by(line_id) %>% 
  dplyr::filter(row_number() == 1)

#Test the significance of the failure rate
test = fisher.test(matrix(c(20,15,122,16), nrow = 2, byrow = TRUE)) #Test accorss differnetiation
test = fisher.test(matrix(c(15,7,115,8), nrow = 2, byrow = TRUE)) #Test accorss lines

#Identify lines that failed to dfferentiate twice
failed_twice = dplyr::filter(diff_fail_count, diff_fail_count > 1)$line_id
dplyr::filter(full_data, line_id %in% failed_twice) %>% arrange(ips_started)

#Look at the lines differentiated in Jan 2015
dplyr::mutate(full_data, days_from = as.numeric(ips_started - as.Date("2015-01-15"))) %>% 
  dplyr::filter(abs(days_from) < 15 ) %>% 
  dplyr::select(line_id, flow_success, ips_started) %>% 
  dplyr::mutate(failed = ifelse(line_id %in% failed_twice, TRUE,FALSE))

#Is there a difference in success rate between feeder and feeder free lines
line_media_assignment = dplyr::select(full_data, line_id, media) %>% group_by(line_id) %>% dplyr::filter(row_number() == 1)
success_rate = dplyr::left_join(diff_fail_count, line_media_assignment, by = "line_id") %>% 
  dplyr::mutate(has_succeeded = diff_attempt_count > diff_fail_count)