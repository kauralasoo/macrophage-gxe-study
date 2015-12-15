library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-chromatin/housekeeping//")
library("ggplot2")

#Load sample neames from disk
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
design_matrix = constructDesignMatrix_ATAC(sample_names)

#Import fragment lengths from disk
fragment_lengths = seqUtils::loadFragmentLengths("processed/SL1344/",sample_names)
normalised_fragment_lengths = dplyr::group_by(fragment_lengths,sample_id) %>%
  dplyr::mutate(frequency = count/sum(count))
fragment_length_plot = dplyr::filter(normalised_fragment_lengths, fragment_length < 350) %>% 
  ggplot(aes(x = fragment_length, y = frequency)) + 
  geom_line() + facet_wrap(~sample_id,ncol = 6)
ggsave("results/ATAC/QC/ATAC_fragment_lengths_by_sample.pdf", fragment_length_plot, width = 8, height = 24)

#Calculate the ratio of short to long fragments
short_fragment_count = dplyr::filter(fragment_lengths, fragment_length < 150) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(short_fragment_count = sum(count))
long_fragment_count = dplyr::filter(fragment_lengths, fragment_length >= 150) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(long_fragment_count = sum(count))
fragment_distribution = dplyr::left_join(short_fragment_count, long_fragment_count, by = "sample_id") %>% 
  dplyr::mutate(length_ratio = short_fragment_count/long_fragment_count) %>%
  dplyr::arrange(-length_ratio)
write.table(fragment_distribution, "macrophage-chromatin/data/SL1344/QC_measures/ATAC_short_long_ratio.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE)


#Fragment lengths before and after adapter trimming
trimmed_lengths = read.table("results/ATAC/bima_A_ATAC.fragment_lengts.txt")
untrimmed_lengths = read.table("results/ATAC/bima_A_ATAC.fragment_lengts.untrimmed.txt")
dat = rbind(dplyr::mutate(untrimmed_lengths, trimmed = "no"), dplyr::mutate(trimmed_lengths, trimmed = "yes")) %>% 
  dplyr::filter(V2 < 300)
fragment_lengths_plot = ggplot(dat,aes(x = V2, y = V1, colour = trimmed)) + geom_point() + geom_line()
ggsave("results/ATAC/QC/trimmed_vs_untrimmed_fragment_lengths.pdf", fragment_lengths_plot, width = 6, height = 6)
