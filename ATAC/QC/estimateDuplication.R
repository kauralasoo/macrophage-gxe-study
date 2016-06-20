library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")

#Import sample names
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

data = seqUtils::loadMarkDuplicates("processed/SL1344", sample_names)
data_filtered = dplyr::select(data, sample_id, PERCENT_DUPLICATION) %>%
  dplyr::rename(percent_duplication = PERCENT_DUPLICATION) %>% 
  dplyr::mutate(percent_duplication = as.numeric(as.character(percent_duplication)))

write.table(data_filtered, "macrophage-chromatin/data/SL1344/QC_measures/ATAC_duplication_fraction.txt", sep = "\t", quote = FALSE, row.names = FALSE)
