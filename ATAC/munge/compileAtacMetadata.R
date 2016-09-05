library("dplyr")
library("ggplot2")
library("tidyr")
library("devtools")
load_all("macrophage-gxe-study/housekeeping/")

#Import compiled line metadata
line_data = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") %>%
  dplyr::filter(!is.na(atac_date))

#Import ATAC sample metadata
atac_meta = read.csv("macrophage-gxe-study/data/chromatin/ATAC/ATAC_sample_metadata.csv", 
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

#Add QC metrics calculated from seq data
assigned = read.table("macrophage-gxe-study/data/chromatin/ATAC/QC_measures/ATAC_assigned_fraction.txt", 
                      header = TRUE, stringsAsFactors = FALSE) %>% 
  dplyr::transmute(sample_id, assigned = Assigned, assigned_frac)
mt_fraction = read.table("macrophage-gxe-study/data/chromatin/ATAC/QC_measures/ATAC_mitochondrial_fraction.txt",
           header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::transmute(sample_id, mt_frac = MT)
dupl_fraction = read.table("macrophage-gxe-study/data/chromatin/ATAC/QC_measures/ATAC_duplication_fraction.txt",
                           header = TRUE, stringsAsFactors = FALSE)
short_long_ratio = read.table("macrophage-gxe-study/data/chromatin/ATAC/QC_measures/ATAC_short_long_ratio.txt",
                           header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::transmute(sample_id, short_long_ratio = length_ratio)
peak_count = read.table("macrophage-gxe-study/data/chromatin/ATAC/QC_measures/macs2_peaks_counts.txt",
                        header = TRUE, stringsAsFactors = FALSE)

atac_metadata = dplyr::left_join(atac_metadata, assigned, by = "sample_id") %>%
  dplyr::left_join(mt_fraction, by = "sample_id") %>%
  dplyr::left_join(dupl_fraction, by = "sample_id") %>%
  dplyr::left_join(short_long_ratio, by = "sample_id") %>%
  dplyr::left_join(peak_count, by = "sample_id")

#Save results to disk
saveRDS(atac_metadata, "macrophage-gxe-study/data/chromatin/ATAC/compiled_atac_metadata.rds")
write.table(atac_metadata, "macrophage-gxe-study/data/chromatin/ATAC/compiled_atac_metadata.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)

