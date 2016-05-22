library("ggplot2")
library("dplyr")
library(ggthemes)

#Compile Pi1 data into a single df
rna_pi1_matrix = read.table("results/SL1344/eQTLs/properties/pi1_matrix_tidy.txt", stringsAsFactors = FALSE) %>% 
  dplyr::mutate(trait = "RNA-seq")
atac_pi1_matrix = read.table("../macrophage-chromatin/results/ATAC/QTLs/properties/fastqtl_pi1_results_tidy.txt", 
                             stringsAsFactors = FALSE) %>%
  dplyr::mutate(trait = "ATAC-seq")
pi1_df = dplyr::bind_rows(rna_pi1_matrix, atac_pi1_matrix) %>% 
  dplyr::mutate(first = factor(first, levels = c("naive","IFNg","SL1344", "IFNg_SL1344"))) %>%
  dplyr::mutate(second = factor(second, levels = rev(c("naive","IFNg","SL1344", "IFNg_SL1344")))) %>%
  dplyr::mutate(trait = factor(trait, levels = c("RNA-seq", "ATAC-seq")))


#Make a heatmap of pi1 values
pi1_plot = ggplot2::ggplot(pi1_df, aes(x = first, y = second, fill = pi1, label = round(pi1,2))) + 
  geom_tile() +
  geom_text() +
  facet_wrap(~trait) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(space = "Lab", low = "#FFFFBF",
                       high = "#E24C36", name = "Pi1") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave("results/SL1344/eQTLs/properties/pi1_heatmap.pdf", plot = pi1_plot, width = 8, height = 4)

  