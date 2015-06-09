library("dplyr")

bioanalyser_dat = read.table("macrophage-gxe-study/data/bioanalyzer/bioanalyser_combined_results.txt", header = TRUE, stringsAsFactors = FALSE)
nanodrop_dat = read.table("macrophage-gxe-study/data/bioanalyzer/matching_nanodrop_data.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE) %>% 
  dplyr::select(sample_id, ng.ul) %>% dplyr::rename(nanodrop = ng.ul)
old_nanodrop_dat = read.table("macrophage-gxe-study/data/RNA_extractions/RNA_extractions_070515.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE) %>%
dplyr::select(sample_id, ng.ul) %>% dplyr::rename(nanodrop_old = ng.ul)

#Merge data:
merged = dplyr::left_join(bioanalyser_dat, old_nanodrop_dat, by = "sample_id") %>% 
  dplyr::left_join(nanodrop_dat, by = "sample_id")