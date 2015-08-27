nanodrop = read.table("macrophage-gxe-study/data/RNA_extractions/RNA_extractions_270815.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
bioanalyzer = read.table("macrophage-gxe-study/data/bioanalyzer/bioanalyzer_samples_2.txt", header = TRUE, stringsAsFactors = FALSE)

combined = dplyr::select(nanodrop, sample_id, ng.ul) %>% dplyr::left_join(bioanalyzer, by = "sample_id")
write.table(combined, "results/SL1344/nano_bio.txt", sep = "\t", quote = FALSE)
