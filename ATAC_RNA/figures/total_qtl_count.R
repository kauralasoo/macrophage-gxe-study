library("purrr")
library("devtools")
library("dplyr")
library("ggplot2")
library("tidyr")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Import QTL counts from each method
caqtl_counts = read.table("results/ATAC/QTLs/properties/fastQTL_vs_rasqual_QTL_counts.txt", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(-extra_qtls) %>%
  tidyr::gather(method, qtl_count, fastqtl_qtl_count, rasqual_qtl_count) %>% 
  dplyr::mutate(method = ifelse(method == "fastqtl_qtl_count", "FastQTL", "RASQUAL")) %>%
  dplyr::mutate(phenotype = "ATAC-seq")
eqtl_counts = read.table("figures/supplementary/rna_fastQTL_vs_rasqual_eQTL_counts.txt", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(-extra_qtls) %>%
  tidyr::gather(method, qtl_count, fastqtl_qtl_count, rasqual_qtl_count) %>% 
  dplyr::mutate(method = ifelse(method == "fastqtl_qtl_count", "FastQTL", "RASQUAL")) %>%
  dplyr::mutate(phenotype = "RNA-seq")

#Merge counts
qtl_counts = dplyr::bind_rows(eqtl_counts, caqtl_counts) %>% 
  dplyr::rename(condition_name = condition) %>%
  dplyr::left_join(housekeeping::figureNames(), by = "condition_name") %>%
  dplyr::mutate(phenotype = factor(phenotype, levels = c("RNA-seq", "ATAC-seq")))

#Make a barplot with the numbers of QTLs
qtl_plot = ggplot(qtl_counts, aes(x = figure_name, y = qtl_count, fill = method)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~phenotype) + 
  theme_light() + 
  ylab("Number of QTLs") + 
  xlab("Condition") +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("#ca0020","#404040"))
ggsave("figures/main_figures/qtl_count_total.pdf", plot = qtl_plot, width = 4, height = 3.5)

