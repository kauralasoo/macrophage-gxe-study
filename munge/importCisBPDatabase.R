library("dplyr")
library("DESeq2")
library("devtools")
library("edgeR")
library("rtracklayer")
load_all("macrophage-chromatin/housekeeping/")
load_all("../seqUtils/")
library("TFBSTools")

#Import TF information form disk
TF_information = readr::read_tsv("~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/TF_Information.txt")
colnames(TF_information)[6] = "gene_id"
unique_motifs = dplyr::select(TF_information, Motif_ID, gene_id, TF_Name) %>% dplyr::filter(Motif_ID != ".")

read.table("~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/pwms_all_motifs/M0082_1.02.txt", header = TRUE)

a = cisBPImportMotif("M0082_1.02","~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/pwms_all_motifs/")


motif_record = unique_motifs[10,]
