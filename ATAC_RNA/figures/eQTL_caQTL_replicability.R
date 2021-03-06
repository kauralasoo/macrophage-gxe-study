library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

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
  calculatePairwiseConcordance(atac_rasqual_50kb, ~p_nominal < 1e-6, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "RASQUAL (50kb)", class = "caQTL")
atac_fastqtl_50kb_concordance = 
  calculatePairwiseConcordance(atac_fastqtl_50kb, ~p_nominal < 1e-6, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (50kb)", class = "caQTL")
atac_fastqtl_100kb_concordance = 
  calculatePairwiseConcordance(atac_fastqtl_100kb, ~p_nominal < 1e-6, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (100kb)", class = "caQTL")

#Calculate RNA concordances
rna_rasqual_concordance = 
  calculatePairwiseConcordance(rna_rasqual_500kb, ~p_nominal < 1e-6, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "RASQUAL (500kb)", class = "eQTL")
rna_fastqtl_100kb_concordance = 
  calculatePairwiseConcordance(rna_fastqtl_100kb, ~p_nominal < 1e-6, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (100kb)", class = "eQTL")
rna_fastqtl_500kb_concordance = 
  calculatePairwiseConcordance(rna_fastqtl_500kb, ~p_nominal < 1e-6, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (500kb)", class = "eQTL")

#Merge all comparisons into a single df
concordance_df = dplyr::bind_rows(rna_rasqual_concordance, rna_fastqtl_100kb_concordance,
                                  rna_fastqtl_500kb_concordance, atac_rasqual_concordance,
                                  atac_fastqtl_50kb_concordance, atac_fastqtl_100kb_concordance)
saveRDS(concordance_df, "results/ATAC_RNA_overlaps/concordance_df.rds")
concordance_df = readRDS("results/ATAC_RNA_overlaps/concordance_df.rds")

#Make a clean concordance df
concordance_df_clean = tidyr::separate(concordance_df, method, into = c("method", "window"), sep = " ") %>% 
  dplyr::mutate(phenotype = ifelse(class == "caQTL", "ATAC-seq", "RNA-seq"))  %>%
  dplyr::mutate(phenotype = factor(phenotype, levels = c("RNA-seq", "ATAC-seq")))

#Make a concordance plot
concordance_filtered = dplyr::filter(concordance_df_clean, window != "(100kb)")
concordance_plot = ggplot(concordance_filtered, aes(x = method, y = concordance, color = method)) + 
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.2)) + 
  facet_wrap(~phenotype) +
  ylab("Lead variant sharing") +
  scale_color_manual(values = c("#ca0020","#404040")) +
  scale_y_continuous(limits = c(0,1)) + 
  theme_light() +
  theme(legend.position = "top", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(colour = "grey10"), strip.background = element_rect(fill = "grey85"))
ggsave("figures/main_figures/concordance_between_conditions_atac_rna.pdf", plot = concordance_plot, width = 4, height = 3.5)

#Make a concordance plot in 100kb cis window
concordance_filtered_100kb = dplyr::filter(concordance_df_clean, window == "(100kb)") %>%
  dplyr::mutate(window = gsub("(", "", window, fixed = TRUE)) %>%
  dplyr::mutate(window = gsub(")", "", window, fixed = TRUE))
concordance_plot_100kb = ggplot(concordance_filtered_100kb, aes(x = window, y = concordance)) + 
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.2)) + 
  facet_wrap(~phenotype) +
  ylab("Lead variant sharing") +
  xlab("Cis-window size") +
  scale_y_continuous(limits = c(0,1)) + 
  theme_light() + 
  theme(strip.text.x = element_text(colour = "grey10"), strip.background = element_rect(fill = "grey85"))
ggsave("figures/supplementary/concordance_between_conditions_atac_rna_100kb_window.pdf", plot = concordance_plot_100kb, width = 3, height = 3.5)



#Calculate Pi1 replicability values
rna_fastqtl_100kb_pi1 = calculatePairwisePi1(rna_fastqtl_100kb, tidy = T) %>%
  dplyr::filter(first != second) %>% 
  dplyr::mutate(window_size = "100kb", class = "RNA-seq")
rna_fastqtl_500kb_pi1 = calculatePairwisePi1(rna_fastqtl_500kb, tidy = T) %>%
  dplyr::filter(first != second) %>% 
  dplyr::mutate(window_size ="500kb", class = "RNA-seq")

atac_fastqtl_100kb_pi1 = calculatePairwisePi1(atac_fastqtl_100kb, tidy = T) %>%
  dplyr::filter(first != second) %>% 
  dplyr::mutate(window_size = "100kb" , class = "ATAC-seq")
atac_fastqtl_50kb_pi1 = calculatePairwisePi1(atac_fastqtl_50kb, tidy = T) %>%
  dplyr::filter(first != second) %>% 
  dplyr::mutate(window_size = "50kb", class = "ATAC-seq")

#Combine into a df and make a plot
pi1_df = dplyr::bind_rows(rna_fastqtl_100kb_pi1, rna_fastqtl_500kb_pi1,
                                  atac_fastqtl_50kb_pi1, atac_fastqtl_100kb_pi1) %>%
  dplyr::mutate(window_size = factor(window_size, levels = c("50kb", "100kb", "500kb")))
pi1_plot = ggplot(pi1_df, aes(x = window_size, y = pi1)) + 
  facet_wrap(~class, scales = "free_x") + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + theme_light() +
  ylab(expression(paste("Pairwise ", pi[1], " statistic"))) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Cis-window size")
ggsave("figures/supplementary/pi1_between_conditions_atac_rna.pdf", pi1_plot, width = 4, height = 3.5)

