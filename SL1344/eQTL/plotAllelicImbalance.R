library("devtools")
library("dplyr")
library("ggplot2")
library("purrr")
library("coloc")
library("readr")
library("optparse")
load_all("../seqUtils/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")
variant_information = importVariantInformation("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

#### ASE data ####
#### NXPH2 ####
#Fetch ASE data from disk
exon_ranges = constructExonRanges("ENSG00000144227", "rs7594476", combined_expression_data$gene_metadata)
sample_meta = dplyr::select(combined_expression_data$sample_metadata, sample_id, condition_name, genotype_id)
ase_data = fetchGeneASEData(exon_ranges, "results/SL1344/combined_ASE_counts.sorted.txt.gz", sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)

#Filter for plotting
plotting_data = dplyr::filter(ase_data, lead_snp_value == 1, feature_snp_value == 1)
ggplot(plotting_data, aes(x = factor(lead_snp_value), y = abs(0.5-ratio), label = sample_id, size = total_count)) + 
  facet_grid(feature_snp_id~condition_name) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = .3, height = .01)) +
  xlab("Feature SNP id") + 
  ylab("Allelic imbalance")

#### SPOPL ####
#Fetch ASE data from disk
exon_ranges = constructExonRanges("ENSG00000144228", "rs7594476", combined_expression_data$gene_metadata)
sample_meta = dplyr::select(combined_expression_data$sample_metadata, sample_id, condition_name, genotype_id)
ase_data = fetchGeneASEData(exon_ranges, "results/SL1344/combined_ASE_counts.sorted.txt.gz", sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)

#Filter for plotting
plotting_data = dplyr::filter(ase_data, lead_snp_value == 1, feature_snp_value == 1)
ggplot(plotting_data, aes(x = factor(lead_snp_value), y = abs(0.5-ratio), label = sample_id)) + 
  facet_grid(feature_snp_id~condition_name) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = .1)) +
  xlab("Feature SNP id") + 
  ylab("Allelic imbalance")


#### SPOPL ####
#Fetch ASE data from disk
exon_ranges = constructExonRanges("ENSG00000185245", "rs4486968", combined_expression_data$gene_metadata)
sample_meta = dplyr::select(combined_expression_data$sample_metadata, sample_id, condition_name, genotype_id)
ase_data = fetchGeneASEData(exon_ranges, "results/SL1344/combined_ASE_counts.sorted.txt.gz", sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)

#Filter for plotting
plotting_data = dplyr::filter(ase_data, lead_snp_value == 1, feature_snp_value == 1)
ggplot(plotting_data, aes(x = factor(lead_snp_value), y = abs(0.5-ratio), label = sample_id)) + 
  facet_grid(feature_snp_id~condition_name) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = .1)) +
  xlab("Feature SNP id") + 
  ylab("Allelic imbalance")

#Make plot
ase_data = fetchGeneASEData(exon_ranges, "results/SL1344/combined_ASE_counts.sorted.txt.gz", sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)
plotting_data = filterASEforPlotting(ase_data) %>% dplyr::filter(lead_snp_value == 1)
new_data_plot = ggplot(plotting_data, aes(x = factor(lead_snp_value), y = abs(0.5-ratio), label = sample_id)) + 
  facet_grid(feature_snp_id~condition_name) +
  geom_boxplot(outlier.shape = NA) + 
  geom_text() +
  geom_jitter(position = position_jitter(width = .1)) +
  xlab("Feature SNP id") + 
  ylab("Allelic imbalance")
ggsave("results/SL1344/SPOPL_AI_new.pdf", new_data_plot, width = 8, height = 10)


new_samples = unique(b$genotype_id)
old_sample_meta = dplyr::filter(sample_meta, !(sample_meta$genotype_id %in% new_samples))

ase_data_old = fetchGeneASEData(exon_ranges, "results/SL1344/combined_ASE_counts.sorted.old.txt.gz", old_sample_meta) %>%
  aseDataAddGenotypes(vcf_file$genotypes)
plotting_data = filterASEforPlotting(ase_data_old) %>% dplyr::filter(lead_snp_value == 1)
old_data_plot = ggplot(plotting_data, aes(x = factor(lead_snp_value), y = abs(0.5-ratio), label = sample_id)) + 
  facet_grid(feature_snp_id~condition_name) +
  geom_boxplot(outlier.shape = NA) +
  geom_text() +
  geom_jitter(position = position_jitter(width = .1)) +
  xlab("Feature SNP id") + 
  ylab("Allelic imbalance")
ggsave("results/SL1344/SPOPL_AI_old.pdf", old_data_plot, width = 8, height = 10)

