library(ggplot2)
library("devtools")
load_all("../seqUtils/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))

#Import rna concentrations
rna_concentrations = readRDS("macrophage-gxe-study/data/covariates/rna_concentrations.rds") %>%
  dplyr::semi_join(combined_expression_data$sample_metadata, by = "sample_id") %>%
  dplyr::select(sample_id, ng_ul)

#Construct a nice df
selected_data = combined_expression_data$sample_metadata %>% 
  dplyr::select(sample_id, condition_name, salmonella_date, donor) %>% 
  dplyr::left_join(rna_concentrations, by = "sample_id") %>% 
  dplyr::group_by(donor) %>% 
  dplyr::mutate(relative_rna = (ng_ul - mean(ng_ul))/sd(ng_ul)) %>%
  dplyr::left_join(friendlyNames(), by = "condition_name")

#Make a plot and save it to disk
relative_rna_plot = ggplot(selected_data, aes(x = friendly_name, y = relative_rna, fill = friendly_name)) + 
  geom_violin() + 
  geom_boxplot(width = 0.3) +
  theme_light() + 
  scale_fill_manual(values = conditionPalette()) + 
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  ylab("Relative RNA amount")
ggsave("results/SL1344/figures/supplementary/relative_rna_amount_per_condition.pdf", 
       relative_rna_plot, width = 4.5, height = 4.5)
ggsave("results/SL1344/figures/supplementary/relative_rna_amount_per_condition.png", 
       relative_rna_plot, width = 4.5, height = 4.5)
