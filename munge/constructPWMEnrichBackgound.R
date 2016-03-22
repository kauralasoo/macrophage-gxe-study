library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("TFBSTools")
library("PWMEnrich")

#Import CisBP motifs
motifs = readRDS("results/ATAC/cisBP_PFMatrixList.rds")
motifs_list = lapply(motifs, Matrix)

cisbp_logn_bg = PWMEnrich::makeBackground(motifs_list, organism = "hg19", algorithm = "human")
saveRDS(cisbp_logn_bg, "results/ATAC/cisBP_PWMEnrich_background.rds")
