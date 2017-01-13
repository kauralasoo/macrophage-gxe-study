library("dplyr")
library("ggplot2")
library("rtracklayer")

#Merge QC data into one table
mt_content = read.table("macrophage-chromatin/data/SL1344/QC_measures/ATAC_mitochondrial_fraction.txt", header = TRUE, stringsAsFactors = FALSE)
sn_ratio = read.table("macrophage-chromatin/data/SL1344/QC_measures/ATAC_assigned_fraction.txt", header = TRUE, stringsAsFactors = FALSE)
duplication_fraction = read.table("macrophage-chromatin/data/SL1344/QC_measures/ATAC_duplication_fraction.txt", header = TRUE, stringsAsFactors = FALSE)
sample_genotype_match = read.table("macrophage-chromatin/data/SL1344/QC_measures/ATAC_sample_genotype_match.txt", header = TRUE, stringsAsFactors = FALSE)
peak_counts = read.table("macrophage-chromatin/data/SL1344/QC_measures/macs2_peaks_counts.txt", header = TRUE, stringsAsFactors = FALSE)
  
#Join all tables together
qc_report = dplyr::left_join(mt_content, sn_ratio, by = "sample_id") %>%
  dplyr::left_join(duplication_fraction, by = "sample_id") %>%
  dplyr::left_join(sample_genotype_match, by = "sample_id") %>%
  dplyr::left_join(peak_counts, by = "sample_id") %>%
  dplyr::select(sample_id, genotype_id, MT, Assigned, assigned_frac, percent_duplication, peak_count)
write.table(qc_report, "macrophage-chromatin/data/SL1344/QC_measures/QC_report.txt", row.names = FALSE, sep ="\t")

#Make histogram of signal-to-noise
sn_plot = ggplot(qc_report, aes(x = assigned_frac)) + 
  geom_histogram(binwidth = 0.05) +
  xlab("Fraction of fragments within peaks") + 
  scale_x_continuous(limits = c(0,1))
ggsave("results/ATAC/QC/ATAC_signal_to_noise.pdf", sn_plot, width = 6, height = 5)

#Make histogram of duplicates
dup_plot = ggplot(qc_report, aes(x = percent_duplication)) + 
  geom_histogram() +
  xlab("Fraction of duplicate fragments") + 
  scale_x_continuous(limits = c(0,0.5))
ggsave("results/ATAC/QC/ATAC_percent_duplication.pdf", dup_plot, width = 6, height = 5)

#Make histogram of assigned fragments
assigned_plot = ggplot(qc_report, aes(x = Assigned/1000000)) + 
geom_histogram(binwidth = 2) +
  xlab("Millions of assigned fragments")
ggsave("results/ATAC/QC/ATAC_assigned_fragments.pdf", assigned_plot, width = 6, height = 5)

#Make histogram of MT fraction
mt_fraction_plot = ggplot(qc_report, aes(x = MT)) + geom_histogram(binwidth = 0.05)
ggsave("results/ATAC/QC/ATAC_mitochondrial_fraction_plot.pdf", mt_fraction_plot, width = 5, height = 5)

#Make histogram of peak counts per sample
peak_count_plot = ggplot(qc_report, aes(x = peak_count)) + geom_histogram(binwidth = 10000)
ggsave("results/ATAC/QC/peak_count_plot.pdf", mt_fraction_plot, width = 5, height = 5)

#### Analyse peak lengths ####
atac_peaks = import.gff3("annotations/chromatin/ATAC_consensus_peaks.gff3")

#Make histogram of peak lengths
peak_length_df = data_frame(length = width(atac_peaks))
plot = ggplot(peak_length_df, aes(x = length)) + 
  geom_histogram() + 
  scale_x_continuous(limits = c(0,1500)) + 
  theme_light() +
  xlab("Peak length (bp)")
ggsave("figures/supplementary/ATAC_peak_length_histogram.png", plot = plot, width = 4, height = 3.5)

#Total peak length
sum(width(atac_peaks))

h3k27ac_peaks = read.table("annotations/H3K27Ac_joint_peaks.bed")
h3k27ac_length = data_frame(peak_length = h3k27ac_peaks$V3-h3k27ac_peaks$V2) %>% dplyr::filter(peak_length < 50000)
