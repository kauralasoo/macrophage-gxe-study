library("readr")
library("devtools")
library("dplyr")
library("tidyr")
library("limma")
library("purrr")
library("ggplot2")
load_all("../seqUtils/")

#Import intron prop data
prop_list = readRDS("results/SL1344/combined_proportions.rds")
prop_list$sample_metadata$condition_name = factor(prop_list$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))

#Perform PCA on the intron proportions
norm_prop_pca = performPCA(prop_list$cqn, prop_list$sample_metadata)
norm_prop_pca$pca_matrix = dplyr::mutate(norm_prop_pca$pca_matrix, protocol = ifelse(rna_auto, "automatic","manual"))

#Make PCA plots
var_exp = norm_prop_pca$var_exp * 100 
xlabel = paste0("PC1 (", round(var_exp[1]), "%)")
ylabel = paste0("PC2 (", round(var_exp[2]), "%)")
xlabel2 = paste0("PC3 (", round(var_exp[3]), "%)")
pca_plot = ggplot(norm_prop_pca$pca_matrix, aes(x = PC1, y = PC3, color = condition_name, shape = protocol)) + 
  geom_point() +
  xlab(xlabel) + 
  ylab(xlabel2) +
  scale_color_manual(values = conditionPalette(), name = "condition") + 
  theme_light() +
  theme(legend.key = element_blank())
ggsave("figures/supplementary/leafcutter_proportions_PCA_PC1_PC3.pdf", plot = pca_plot, width = 5.5, height = 4)

pca_plot2 = ggplot(norm_prop_pca$pca_matrix, aes(x = PC1, y = PC2, color = condition_name, shape = protocol)) + 
  geom_point() +
  xlab(xlabel) + 
  ylab(ylabel) +
  scale_color_manual(values = conditionPalette(), name = "condition") + 
  theme_light() +
  theme(legend.key = element_blank())
ggsave("figures/supplementary/leafcutter_proportions_PCA_PC1_PC2.pdf", plot = pca_plot2, width = 5.5, height = 4)


