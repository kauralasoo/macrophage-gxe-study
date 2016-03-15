library("PWMEnrich")
library("dplyr")

#Import ATAC data and clusters
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))
final_clusters = readRDS("results/ATAC/DA/peak_clusters.rds")


#Focus on Drosophila example
library(PWMEnrich.Dmelanogaster.background)
library("PWMEnrich.Hsapiens.background")
data(PWMLogn.dm3.MotifDb.Dmel)

#Calculate enrichments
sequence = readDNAStringSet(system.file(package="PWMEnrich", dir="extdata", file="stripe2.fa"))
res = motifEnrichment(sequence, PWMLogn.dm3.MotifDb.Dmel)
report = sequenceReport(res, 1)

#Look at binding sites
ids = c("bcd_FlyReg_FBgn0000166",
        "gt_FlyReg_FBgn0001150",
        "Kr")
sel.pwms = PWMLogn.dm3.MotifDb.Dmel$pwms[ids]
names(sel.pwms) = c("bcd", "gt", "Kr")

scores = motifScores(sequence, sel.pwms, raw.scores=TRUE)
plotMotifScores(scores, cols=c("green", "red", "blue"))

#Look for enrichmen in multiple sequences
# load the pre-compiled lognormal background
data(PWMLogn.dm3.MotifDb.Dmel)
sequences = readDNAStringSet(system.file(package="PWMEnrich",
                                         dir="extdata", file="tinman-early-top20.fa"))
res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
report = groupReport(res)

#Import ATAC peak sequences from disk
sequences = readDNAStringSet("annotations/ATAC_Seq_joint_peaks.fasta")
peak_ids = strsplit(names(sequences), "=") %>% lapply(function(x) x[2]) %>% unlist()
names(sequences) = peak_ids

#Human data
data(PWMLogn.hg19.MotifDb.Hsap)
res = motifEnrichment(DNAString("TGCATCAAGTGTGTAGTGCAAGTGAGTGATGAGTAGAAGTTGAGTGAGGTAGATGC"), PWMLogn.hg19.MotifDb.Hsap)
groupReport(res)





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
