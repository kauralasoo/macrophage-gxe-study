library("devtools")
library("dplyr")
library("lme4")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")

