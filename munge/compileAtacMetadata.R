library("dplyr")
library("ggplot2")
library("tidyr")
library("devtools")
load_all("macrophage-chromatin/housekeeping/")

#Import compiled line metadata
line_data = readRDS("../macrophage-gxe-study/macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") %>%
  dplyr::filter(!is.na(atac_date))

#Import ATAC sample metadata
atac_meta = read.csv("macrophage-chromatin/data/SL1344/ATAC_sample_metadata.csv", 
         stringsAsFactors = FALSE, na.strings = "")
atac_design = constructDesignMatrix_ATAC(atac_meta$sample_id)

#Join'em together
atac_metadata = dplyr::left_join(atac_design, atac_meta, by = "sample_id") %>%
  dplyr::mutate(stimulation_date = as.Date(stimulation_date, "%d/%m/%Y")) %>%
  dplyr::mutate(submission_date = as.Date(submission_date, "%d/%m/%Y")) %>%
  dplyr::left_join(line_data, by = "donor") %>%
  dplyr::select(-comment.x, -comment.y)

#Remove flow purity for samples were flow was done more than 2 weeks after RNA extraction
atac_metadata = dplyr::mutate(atac_metadata, 
            max_purity_filtered = ifelse(as.numeric(flow_date - stimulation_date) > 14, NA, max_purity),
            mean_purity_filtered = ifelse(as.numeric(flow_date - stimulation_date) > 14,NA,mean_purity))

#Save results to disk
saveRDS(line_data, "macrophage-chromatin/data/SL1344/compiled_atac_metadata.rds")
write.table(line_data, "macrophage-chromatin/data/SL1344/compiled_atac_metadata.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)

