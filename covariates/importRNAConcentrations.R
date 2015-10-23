library("dplyr")
library("ggplot2")
library("tidyr")

rna_data = tbl_df(read.table("macrophage-gxe-study/data/RNA_extractions/RNA_extractions_070515.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE))

#Select and process relevant columns
rna_df = dplyr::rename(rna_data, ng_ul = ng.ul, ratio_260_280 = X260.280, ratio_260_230 = X260.230) %>% #Rename columns
  dplyr::filter(!is.na(sample_id)) %>% #Keep non-NA samples
  dplyr::select(sample_id, Date, ng_ul, ratio_260_280, ratio_260_230) %>% #Select only relevant columns
  dplyr::mutate(Date = as.Date(Date, "%m/%d/%Y"))

#Deal with know artifacts in the data ruyv_C
rna_sample_stats = dplyr::mutate(rna_df, ng_ul = ifelse(grepl("eofe_",sample_id), 2*ng_ul, ng_ul)) %>% #1 well per condition
  dplyr::mutate(ng_ul = ifelse(grepl("fpdl_",sample_id), 2*ng_ul, ng_ul)) %>% #1 well per condition
  dplyr::mutate(ng_ul = ifelse(sample_id == "gomv_D", (21.5/35*ng_ul), ng_ul)) %>% #gomv_D had 21.5 ul of sample instead of 35 ul
  dplyr::mutate(ng_ul = ifelse(sample_id %in% c("ieki_A", "ieki_B"), 2*ng_ul, ng_ul)) %>% #ieki A,B only one well per condition
  dplyr::mutate(ng_ul = ifelse(sample_id %in% c("ougl_A_2","ougl_B_2","ougl_C_2","ougl_D_2"),2*ng_ul,ng_ul)) %>% #ougl_X_2 one well per condition
  dplyr::mutate(ng_ul = ifelse(sample_id == "ruyv_C",2*ng_ul,ng_ul)) %>% #ougl_X_2 one well per condition 
  dplyr::mutate(ng_ul = ifelse(sample_id %in% c("pelm_A","pelm_B","pelm_C","pelm_D"),2*ng_ul,ng_ul)) #pelm_X one well per condition

#Separate sample_id into fields
rna_sample_stats = tidyr::separate(rna_sample_stats, sample_id, into = c("donor","condition","replicate"), sep = "_", extra = "drop", remove = FALSE) %>%
  mutate(replicate = ifelse(is.na(replicate), 1, replicate)) %>% #if no replicate in file name then set it to 1
  mutate(replicate = as.integer(replicate))
  
#Calculate mean concentration per line
rna_line_stats = dplyr::filter(rna_sample_stats, sample_id != "aipt_A") %>%
  dplyr::group_by(donor, replicate) %>%
  dplyr::summarise(ng_ul_mean = mean(ng_ul), extraction_date = Date[1]) %>%
  dplyr::ungroup()

#Combine line-level and sample-level stats
rna_stats = dplyr::left_join(rna_sample_stats, rna_line_stats, by = c("donor", "replicate"))
saveRDS(rna_stats, "macrophage-gxe-study/data/covariates/rna_concentrations.rds")
write.table(rna_stats, "macrophage-gxe-study/data/covariates/rna_concentrations.txt", sep = "\t", quote = FALSE, row.names = FALSE)

