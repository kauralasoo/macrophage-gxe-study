library("flowCore")
library("flowViz")
library("tidyr")
library("plyr")
library("dplyr")
library("openCyto")
library("data.table")
library("ggplot2")

#Import flowSet from disk
flowset = readRDS("flow/combined_flowSet.rds")
metadata = pData(flowset) %>% tbl_df()

#Read the gating information
gating = read.table("macrophage-gxe-study/data/flow/gating_parameters_2.csv", sep = ",", header = TRUE)
gt = gatingTemplate("macrophage-gxe-study/data/flow/gating_parameters_2.csv")

#Transfrom stained channels
transFuncts = estimateLogicle(flowset[[1]], channels = c("PE.A", "APC.A", "Pacific.Blue.A"))
flowset_trans = transform(flowset, transFuncts)

#Perform gating
gs = GatingSet(flowset_trans)
gating(gt, gs)

#Plot some examples
plotGate(gs[["fpdl_291014_CD14+CD16+CD206"]], path = 1, default.y = "SSC.A")
plotGate(gs[["fpdl_291014_isotype"]], path = 1, default.y = "SSC.A")

#Extract population statistics for QC purposes
stat = t(getPopStats(gs))

#Extract data from different gated populations
ungated_data = getData(gs, "root")
nondebris_data = getData(gs, "nonDebris")
gated_data = getData(gs, "cells")

#Convert flowSet into a gigantic data frame
gated_list = as(gated_data, "list")
df_list = lapply(gated_list, function(flowframe){tbl_df(as.data.frame(exprs(flowframe)))})
gated_df = ldply(df_list,.id = "name") %>% tbl_df() %>% dplyr::mutate(name = as.character(name))

#Prepare data for plotting
selected_df = dplyr::select(gated_df, name, APC.A, PE.A, Pacific.Blue.A) %>%
  tidyr::gather(channel, intensity, APC.A:Pacific.Blue.A) %>% 
  dplyr::left_join(metadata, by = "name")

#Make density plots for all channels
a = dplyr::filter(selected_df, sample == "aipt_160614")
plot = ggplot(selected_df, aes(x = intensity, color = staining, fill = staining, alpha = 0.5)) + geom_density() + facet_grid(sample~channel) 
ggsave("results/flow_density.pdf", plot = plot, width = 8, height = 49)




