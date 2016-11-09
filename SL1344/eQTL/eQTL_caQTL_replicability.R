library("devtools")
library("dplyr")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("../seqUtils/")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import eQTL p-values
rna_fastqtl_500kb = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues_500kb.rds")
rna_fastqtl_100kb = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues_100kb.rds")
rna_rasqual_500kb = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")

#Import caQTL p-valurs
atac_fastqtl_100kb = readRDS("results/ATAC/QTLs/fastqtl_min_pvalues_100kb.rds")
atac_fastqtl_50kb = readRDS("results/ATAC/QTLs/fastqtl_min_pvalues_50kb.rds")
atac_rasqual_50kb = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")

#Calculate chromatin concordances
atac_rasqual_concordance = 
  calculatePairwiseConcordance(atac_rasqual_50kb, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "RASQUAL (50kb)", class = "caQTL")
atac_fastqtl_50kb_concordance = 
  calculatePairwiseConcordance(atac_fastqtl_50kb, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (50kb)", class = "caQTL")
atac_fastqtl_100kb_concordance = 
  calculatePairwiseConcordance(atac_fastqtl_100kb, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (100kb)", class = "caQTL")

#Calculate RNA concordances
rna_rasqual_concordance = 
  calculatePairwiseConcordance(rna_rasqual_500kb, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "RASQUAL (500kb)", class = "eQTL")
rna_fastqtl_100kb_concordance = 
  calculatePairwiseConcordance(rna_fastqtl_100kb, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (100kb)", class = "eQTL")
rna_fastqtl_500kb_concordance = 
  calculatePairwiseConcordance(rna_fastqtl_500kb, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (500kb)", class = "eQTL")

#Merge all comparisons into a single df
concordance_df = dplyr::bind_rows(rna_rasqual_concordance, rna_fastqtl_100kb_concordance,
                                  rna_fastqtl_500kb_concordance, atac_rasqual_concordance,
                                  atac_fastqtl_50kb_concordance, atac_fastqtl_100kb_concordance)
saveRDS(concordance_df, "results/ATAC_RNA_overlaps/concordance_df.rds")
concordance_df = readRDS("results/ATAC_RNA_overlaps/concordance_df.rds")

#Make a concordance plot
concordance_df_clean = dplyr::mutate(concordance_df, 
            class = factor(class, levels = c("eQTL","caQTL")))
concordance_plot = ggplot(concordance_df_clean, aes(x = method, y = concordance)) + 
                                       facet_wrap(~class, scales = "free_x") + 
                                       geom_boxplot() + geom_point() + theme_light() +
                                       ylab("Lead variant concordance") +
                                       scale_y_continuous(limits = c(0,1)) +
                                       theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
ggsave("figures/main_figures/concordance_between_conditions_atac_rna.pdf", concordance_plot, width = 6, height = 4)

#Calculate Pi1 replicability values
rna_fastqtl_100kb_pi1 = calculatePairwisePi1(rna_fastqtl_100kb, tidy = T) %>%
  dplyr::filter(first != second) %>% 
  dplyr::mutate(method = "FastQTL (100kb)", class = "eQTL")
rna_fastqtl_500kb_pi1 = calculatePairwisePi1(rna_fastqtl_500kb, tidy = T) %>%
  dplyr::filter(first != second) %>% 
  dplyr::mutate(method = "FastQTL (500kb)", class = "eQTL")

atac_fastqtl_100kb_pi1 = calculatePairwisePi1(atac_fastqtl_100kb, tidy = T) %>%
  dplyr::filter(first != second) %>% 
  dplyr::mutate(method = "FastQTL (100kb)", class = "caQTL")
atac_fastqtl_50kb_pi1 = calculatePairwisePi1(atac_fastqtl_50kb, tidy = T) %>%
  dplyr::filter(first != second) %>% 
  dplyr::mutate(method = "FastQTL (50kb)", class = "caQTL")

#Combine into a df and make a plot
pi1_df = dplyr::bind_rows(rna_fastqtl_100kb_pi1, rna_fastqtl_500kb_pi1,
                                  atac_fastqtl_50kb_pi1, atac_fastqtl_100kb_pi1) %>%
  dplyr::mutate(class = factor(class, levels = c("eQTL","caQTL")))
pi1_plot = ggplot(pi1_df, aes(x = method, y = pi1)) + 
  facet_wrap(~class, scales = "free_x") + 
  geom_boxplot() + geom_jitter(width = 0.1) + theme_light() +
  ylab(expression(paste("Pairwise ", pi[1], " statistic"))) +
  scale_y_continuous(limits = c(0,1)) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
ggsave("figures/main_figures/pi1_between_conditions_atac_rna.pdf", pi1_plot, width = 6, height = 4)

