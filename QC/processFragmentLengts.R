fragment_table = read.table("processed/SL1344/bima_A_ATAC/bima_A_ATAC.fragments_lengts.txt")
fragment_table = dplyr::filter(fragment_table, length < 500)
colnames(fragment_table) = c("count", "length")

fragment_plot = ggplot(fragment_table, aes(x = length, y = count)) + geom_point() + geom_line()
ggsave("results/ATAC/QC/ATAC_fragment_length_plot.pdf", fragment_plot, width = 8, height  =5)
