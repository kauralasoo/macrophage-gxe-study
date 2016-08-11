library("devtools")
library("dplyr")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")

#Import read count data
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Extract PEER factors for different conditions
condition_names = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
rna_conditions = lapply(condition_names, extractConditionFromExpressionList, combined_expression_data)
sample_metadata_list = lapply(rna_conditions, function(x){x$sample_metadata})

#Focus on the relationship between purity and first PEER factor
purity_df = dplyr::filter(sample_metadata_list$naive, mean_purity_filtered > 0.90)
peer1_plot = ggplot(purity_df, aes(x = mean_purity_filtered, y = PEER_factor_1)) + 
  geom_point() +
  xlab("Mean purity") + 
  ylab("PEER factor 1") +
  theme_light()
ggsave("figures/supplementary/varComp_purity_vs_PEER_factor_1.pdf", peer1_plot, width = 4, height = 4)
cor(purity_df$mean_purity_filtered, purity_df$PEER_factor_1, use = "complete.obs")

#RNA concetration vs PEER 1
peer1_plot_2 = ggplot(purity_df, aes(x = ng_ul_mean, y = PEER_factor_1)) + 
  geom_point() +
  xlab("RNA concentration (ng/ul)") + 
  ylab("PEER factor 1") + 
  theme_light()
ggsave("figures/supplementary/varComp_rna_concetration_vs_PEER_factor_1.pdf", peer1_plot_2, width = 4, height = 4)
cor.test(purity_df$ng_ul_mean, purity_df$PEER_factor_1, use = "complete.obs")

#Focus on the relationship between sequencing chemistry and first PEER factor
peer2_plot = ggplot(sample_metadata_list$naive, aes(x = rna_auto, y = PEER_factor_2)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = .1)) + 
  xlab("Automatic RNA-seq") + 
  ylab("PEER factor 2") +
  theme_light()
ggsave("figures/supplementary/varComp_autmatic_vs_PEER_factor_2.pdf", peer2_plot, width = 4, height = 4)
cor.test(as.numeric(sample_metadata_list$naive$rna_auto) , sample_metadata_list$naive$PEER_factor_2, use = "complete.obs")

#### Look at the variation explained by each PEER factor ####
prec_naive = read.table("results/SL1344/PEER/naive_10/precision.txt", sep = ",")
prec_ifng = read.table("results/SL1344/PEER/IFNg_10/precision.txt", sep = ",")
prec_sl1344 = read.table("results/SL1344/PEER/SL1344_10//precision.txt", sep = ",")
prec_ifng_sl1344 = read.table("results/SL1344/PEER/IFNg_SL1344_10//precision.txt", sep = ",")

#Combine it all together
precision_df = 1/cbind(prec_naive, prec_ifng, prec_sl1344, prec_ifng_sl1344)[-1,]
colnames(precision_df) = c("naive", "IFNg", "SL1344", "IFNg_SL1344")
precision_df = dplyr::mutate(precision_df, factor = paste("PEER_factor_", c(1:10), sep = "")) %>% 
  dplyr::mutate(factor = factor(factor, factor)) %>%
  tidyr::gather(condition, sd, naive:IFNg_SL1344) %>%
  dplyr::group_by(condition) %>% 
  dplyr::mutate(var_exp = sd/sum(sd))

#Make a plot
precision_plot = ggplot(precision_df, aes(x = factor, y = var_exp, color = condition, group = condition)) +
  geom_point() + 
  geom_line() +
  ylab("Fraction variance explained") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 25))
precision_plot
ggsave("figures/supplementary/varComp_PEER_explained.pdf", precision_plot, width = 6, height = 5)



