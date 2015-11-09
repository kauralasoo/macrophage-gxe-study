library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")

sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

data = seqUtils::loadFeaturCountsSummary("processed/SL1344", sample_names, counts_suffix = ".consensus_peaks.counts.txt.summary")
data_filtered = dplyr::select(data, sample_id, Assigned, Unassigned_Ambiguity, Unassigned_NoFeatures) %>%
  dplyr::mutate(assigned_frac = Assigned/(Assigned + Unassigned_NoFeatures)) %>%
  dplyr::mutate(ambiguity_frac = Unassigned_Ambiguity/(Assigned + Unassigned_Ambiguity))

#Save to disk
write.table(data_filtered, "macrophage-chromatin/data/SL1344/QC_measures/ATAC_assigned_fraction.txt", sep = "\t", quote = FALSE, row.names = FALSE)

