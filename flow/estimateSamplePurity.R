library("flowCore")
library("flowViz")
library("tidyr")
library("plyr")
library("dplyr")
library("openCyto")
library("ggplot2")
library("mclust")

#Import flowSet from disk
flowset = readRDS("flow/combined_flowSet.rds")
metadata = pData(flowset) %>% tbl_df()

#Filter out QC-fail samples
qc_fail_samples = c("iasn_160614","liun_160614","piun_271114","piun-SN_271114","qanu_160614")
filtered_metadata = dplyr::filter(metadata, !(sample %in% qc_fail_samples))
flowset = flowset[filtered_metadata$name]

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
plotGate(gs[["ffdk_200514_CD14+CD16+CD206"]], path = 1, default.y = "SSC.A")
plotGate(gs[["ougl_200514_CD14+CD16+CD206"]],path = 1, default.y = "SSC.A")

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
plot = ggplot(selected_df, aes(x = intensity, color = staining, fill = staining, alpha = 0.5)) + geom_density() + facet_grid(sample~channel) 
ggsave("results/flow_density.pdf", plot = plot, width = 8, height = 49)

#Estimate the purity of all samples
purity_data = ldply(as.list(sample_names), estimateSamplePurity, filtered_metadata, gated_data)
purity_data = tidyr::separate(purity_data, sample, c("donor", "date"), sep ="_", remove = FALSE)
write.table(purity_data, "results/flow/purity_estimates.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#### Finally, we want to know if the purity values are lower for those lines that cluster separately on the RNA-Seq PCA plot. ####

#Keep only the lines that also have RNA-Seq data
rna_design = readRDS("results/SL1344/design_matrix.rds")
cluster_samples = c("iasn","huls","debk", "ougl", "gomv", "ffdp", "peop")
outliers = c("golb_111114","fpdj_200514")

#Calculate and filter max purity
max_purity = dplyr::group_by(purity_data, sample) %>% summarize(max_purity = max(purity), donor = donor[1])
filtered_purity = semi_join(max_purity, rna_design, by = "donor")
filtered_purity = dplyr::mutate(filtered_purity, cluster = ifelse(donor %in% cluster_samples, "yes", "no")) %>%
  dplyr::filter(!(sample %in% outliers))

#Make plot
plot = ggplot(filtered_purity, aes(x = cluster, y = max_purity, label = donor)) + geom_boxplot() + geom_point() + geom_text()
ggsave("results/flow/purity_difference.pdf", plot = plot)


