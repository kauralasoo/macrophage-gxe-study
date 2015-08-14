library("dplyr")
library("ggplot2")

loadChrCounts <- function(sample_dir, sample_names, counts_suffix = ".chr_counts"){
  #Import read counts per chromosome for each sample
  matrix = c()
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    print(path)
    table = read.table(path, header = FALSE, stringsAsFactors = FALSE)
    colnames(table) = c(sample_names[i], "chr_name")
    if (i == 1){
      matrix = table[,c(2,1)]
      print(head(matrix))
    }
    else{
      matrix = dplyr::full_join(matrix, table, by = "chr_name")
      print(head(matrix))
    }
  }
  return(matrix)
}

#Load sample neames from disk
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

#Load data from disk
data = loadChrCounts("processed/SL1344", sample_names, counts_suffix = ".chr_counts")
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
