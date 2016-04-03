library(Biostrings)
library("TFBSTools")

#Perform analysis with TFBSTools
cisbp_pwm_list = readRDS("results/ATAC/cisBP_PWMatrixList.rds")

#Import ATAC peak sequences from disk
sequences = readDNAStringSet("annotations/ATAC_consensus_peaks.fasta")
peak_ids = strsplit(names(sequences), "=") %>% lapply(function(x) x[2]) %>% unlist()
names(sequences) = peak_ids

#Search for matches
sitesetList = searchSeq(cisbp_pwm_list, a[1:2], min.score="80%", strand="*")
b = as.data.frame(sitesetList)


