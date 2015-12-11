trimmed_lengths = read.table("results/ATAC/bima_A_ATAC.fragment_lengts.txt")
untrimmed_lengths = read.table("results/ATAC/bima_A_ATAC.fragment_lengts.untrimmed.txt")
dat = rbind(dplyr::mutate(untrimmed_lengths, trimmed = "no"), dplyr::mutate(trimmed_lengths, trimmed = "yes")) %>% 
  dplyr::filter(V2 < 300)
fragment_lengths_plot = ggplot(dat,aes(x = V2, y = V1, colour = trimmed)) + geom_point() + geom_line()
ggsave("results/ATAC/QC/trimmed_vs_untrimmed_fragment_lengths.pdf", fragment_lengths_plot, width = 6, height = 6)
