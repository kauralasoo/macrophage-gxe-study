library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")
library("readr")
library("limma")
library("edgeR")
library("Biobase")
library("MFuzz")
library("ggplot2")

#Import data
counts = readRDS("results/Ivashkiv/STAT1_combined_counts.rds")
stat1_names = c("STAT1_rep1_A","STAT1_rep1_B","STAT1_rep1_C","STAT1_rep1_D")
stat1_design = constructDesignMatrix_Ivashkiv(stat1_names)
stat1_peaks = import.gff3("annotations/STAT1_joint_peaks.gff3")

#Construct count matrix
rownames(counts) = counts$gene_id
count_matrix = counts[,stat1_design$sample_id]
length_df = dplyr::select(counts, gene_id, length)
stat1_tpm = log(calculateTPM(count_matrix, length_df) + 0.01,2) #Take the log
stat1_tpm = stat1_tpm[!(apply(stat1_tpm, 1, min) < 0),] #Filter out some akward peaks

#Pick the most variable peaks
stat1_variable_peaks = stat1_tpm[apply(stat1_tpm, 1, sd) > 0.5,]
stat1_set = ExpressionSet(as.matrix(stat1_variable_peaks))
stat1_set_std = standardise(stat1_set)

#Cluster the expression data
clusters = mfuzz(stat1_set_std, c = 7, m = 2.5, iter = 100)
cluster_cores = acore(cqn_set_std, clusters, min.acore = 0.7)
names(cluster_cores) = c(1:length(cluster_cores))

pheatmap(clusters$centers)
