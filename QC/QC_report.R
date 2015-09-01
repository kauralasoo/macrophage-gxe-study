#Merge QC data into one table
mt_content = read.table("results/ATAC/QC/ATAC_mitochondrial_fraction.txt", header = TRUE, stringsAsFactors = FALSE)
sn_ratio = read.table("results/ATAC/QC/ATAC_assigned_fraction.txt", header = TRUE, stringsAsFactors = FALSE)
duplication_fraction = read.table("results/ATAC/QC/ATAC_duplication_fraction.txt", header = TRUE, stringsAsFactors = FALSE)
sample_genotype_match = read.table("results/ATAC/QC/ATAC_sample_genotype_match.txt", header = TRUE, stringsAsFactors = FALSE)

#Join all tables together
qc_report = dplyr::left_join(mt_content, sn_ratio, by = "sample_id") %>%
  dplyr::left_join(duplication_fraction, by = "sample_id") %>%
  dplyr::left_join(sample_genotype_match, by = "sample_id") %>%
  dplyr::select(sample_id, genotype_id, MT, Assigned, assigned_frac, percent_duplication)

#Make histogram of signal-to-noise
sn_plot = ggplot(qc_report, aes(x = assigned_frac)) + 
  geom_histogram(binwidth = 0.02) +
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
geom_histogram(binwidth = 1) +
  xlab("Millions of assigned fragments")
ggsave("results/ATAC/QC/ATAC_assigned_fragments.pdf", assigned_plot, width = 6, height = 5)

#### Analyse peak lengths ####
atac_peaks = read.table("annotations/ATAC_Seq_joint_peaks.gff3")
peak_length = data.frame(peak_length = atac_peaks$V3-atac_peaks$V2)
#Total peak length
sum(peak_length$peak_length)

h3k27ac_peaks = read.table("annotations/H3K27Ac_joint_peaks.bed")
h3k27ac_length = data_frame(peak_length = h3k27ac_peaks$V3-h3k27ac_peaks$V2) %>% dplyr::filter(peak_length < 50000)
