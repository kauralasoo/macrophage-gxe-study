library("rtracklayer")
library("devtools")
library("dplyr")
load_all("../macrophage-gxe-study/macrophage-gxe-study/seqUtils/")

loadNarrowPeaks <- function(sample_dir, sample_names, peaks_suffix = "_peaks.narrowPeak"){
  #Import read counts per chromosome for each sample
  result = list()
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], peaks_suffix, sep = ""))
    print(path)
    peaks = rtracklayer::import(path, format = "BED", extraCols = extraCols_narrowPeak)
    result[[sample_names[i]]] = peaks
  }
  return(result)
}

filterOverlaps <- function(peak_list, minOverlapCount = 3){
  #For each set of peaks, find the ones that are present in at least minOverlapCount samples
  result = list()
  for (i in 1:length(peak_list)){
    query_peaks = peak_list[[i]]
    query_hits = c()
    for (j in 1:length(peak_list)){
      target_peaks = peak_list[[j]]
      overlaps = findOverlaps(query_peaks, target_peaks)
      unique_hits = unique(queryHits(overlaps))
      query_hits = c(query_hits, unique_hits)
    }
    query_hits = table(query_hits)
    query_hits = query_hits[query_hits >= minOverlapCount]
    query_filtered = query_peaks[as.numeric(rownames(query_hits))]
    result[i] = query_filtered
  }
  return(result)
}

#Find peaks in condition A
cond_A_names = c("bima_A_ATAC","vass_A_ATAC","cicb_A_ATAC","eofe_A_ATAC")
cond_A_peaks = loadNarrowPeaks("processed/SL1344/", cond_A_names)
cond_A_filtered = filterOverlaps(cond_A_peaks, 3)
cond_A_union = listUnion(cond_A_filtered)

#Find peaks in condition B
cond_B_names = c("bima_B_ATAC","vass_B_ATAC","cicb_B_ATAC","eofe_B_ATAC")
cond_B_peaks = loadNarrowPeaks("processed/SL1344/", cond_B_names)
cond_B_filtered = filterOverlaps(cond_B_peaks, 3)
cond_B_union = listUnion(cond_B_filtered)

#Find peaks in condition C
cond_C_names = c("bima_C_ATAC","vass_C_ATAC","cicb_C_ATAC","eofe_C_ATAC")
cond_C_peaks = loadNarrowPeaks("processed/SL1344/", cond_C_names)
cond_C_filtered = filterOverlaps(cond_C_peaks, 3)
cond_C_union = listUnion(cond_C_filtered)

#Find peaks in condition D
cond_D_names = c("bima_D_ATAC","vass_D_ATAC","cicb_D_ATAC","eofe_D_ATAC")
cond_D_peaks = loadNarrowPeaks("processed/SL1344/", cond_D_names)
cond_D_filtered = filterOverlaps(cond_D_peaks, 3)
cond_D_union = listUnion(cond_D_filtered)

#Join all peaks together
all_peaks = listUnion(list(cond_A_union, cond_B_union, cond_C_union, cond_D_union))
metadata = data_frame(type = "exon", gene_id = paste("ATAC_peak_", c(1:length(all_peaks)), sep = ""))
elementMetadata(all_peaks) = metadata

#Export as a gff file
rtracklayer::export.gff3(all_peaks, "annotations/ATAC_Seq_joint_peaks.gff3")
rtracklayer::export.bed(all_peaks, "annotations/ATAC_Seq_joint_peaks.bed")
