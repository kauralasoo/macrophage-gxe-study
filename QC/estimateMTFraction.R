library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")

#Load sample neames from disk
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

#Load data from disk
data = seqUtils::loadChrCounts("processed/SL1344", sample_names, counts_suffix = ".chr_counts")
rownames(data) = data$chr_name
data_matrix = data[,-1]

#Estimate the fraction of reads from mitochondria and unmapped reads (possibly salmonella?)
mt_fraction = as.data.frame(t(data_matrix[c("MT","*"),]/colSums(data_matrix, na.rm = TRUE)))
colnames(mt_fraction) = c("MT", "unmapped")
mt_data = dplyr::mutate(mt_fraction, sample_id = rownames(mt_fraction)) %>%
  dplyr::select(sample_id, everything())
write.table(mt_data, "results/ATAC/QC/ATAC_mitochondrial_fraction.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Make histogram
mt_fraction_plot = ggplot(mt_data, aes(x = MT)) + geom_histogram(binwidth = 0.05)
ggsave("results/ATAC/QC/ATAC_mitochondrial_fraction_plot.pdf", mt_fraction_plot, width = 5, height = 5)


