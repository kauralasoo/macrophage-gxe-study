library("PWMEnrich")

#Import ATAC data and clusters
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))
final_clusters = readRDS("results/ATAC/DA/peak_clusters.rds")















#Play around with TFBS_DB
#Perform motif enrichment
db <- CisBP.extdata("Homo_sapiens")
tfs <- tfbs.createFromCisBP(db)

#Agnes clustering
tfs1 <- tfbs.clusterMotifs(tfs, method="agnes", group.k=100, pdf.heatmap="agnes.hm.pdf")
#AP clustering
tfs2 <- tfbs.clusterMotifs(tfs, method="apcluster", pdf.heatmap= "apcluster.hm.pdf")

# Draw motif logos with one group of TF per page
tfbs.drawLogosForClusters(tfs1, file.pdf="agnes.logos.pdf");
tfbs.drawLogosForClusters(tfs2, file.pdf="apcluster.logos.pdf")
