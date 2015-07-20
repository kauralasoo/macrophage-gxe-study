library(dplyr)

id_map = read.table("fastq/id_map_2.txt", stringsAsFactors = FALSE)
colnames(id_map) = c("sequencescape", "sample_id")

sample_names = read.table("fastq/SL1344_names_3.txt", stringsAsFactors = FALSE, comment.char = "")
colnames(sample_names) = c("sequencescape", "lanelets")

new_sample_names = dplyr::left_join(id_map, sample_names, by = "sequencescape") %>% dplyr::select(-sequencescape)
write.table(new_sample_names, "macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
