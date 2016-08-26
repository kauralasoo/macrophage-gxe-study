library("devtools")
library("dplyr")
library("purrr")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")

#Import read count data
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Extract PEER factors for different conditions
condition_names = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
rna_conditions = lapply(condition_names, extractConditionFromExpressionList, combined_expression_data)
sample_metadata_list = lapply(rna_conditions, function(x){x$sample_metadata})

#Look at the distribution of total variance in each condition
cqn_list = purrr::map(rna_conditions,~.$cqn)
filtered_cqn_list = purrr::map(cqn_list, ~.[rowMeans(.) > 2,])
standrard_variance = purrr::map(filtered_cqn_list, function(x){data_frame(var = apply(x,1,sd)/rowMeans(x))}) %>%
  purrr::map_df(identity, .id = "condition_name") %>%
  dplyr::filter(var < 0.5) %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
#Make a plot
variance_per_condition = ggplot(standrard_variance, aes(x = var, color = condition_name)) + 
  geom_density() + 
  theme_light() +
  scale_color_manual(values = conditionPalette()) + 
  xlab("Coefficient of variation")
ggsave("figures/supplementary/varComp_CV_distribution.pdf", plot = variance_per_condition, width = 4.5, height = 3.5)

#Focus on the relationship between purity and first PEER factor
purity_df = dplyr::filter(sample_metadata_list$naive, mean_purity_filtered > 0.90)
purity_test = cor.test(purity_df$mean_purity_filtered, purity_df$PEER_factor_1, use = "complete.obs")
purity_text = paste0("r=",round(purity_test$estimate, 2),", " ,"p=", signif(purity_test$p.value,3))
peer1_plot = ggplot(purity_df, aes(x = mean_purity_filtered, y = PEER_factor_1)) + 
  geom_point() +
  xlab("Mean purity") + 
  ylab("PEER factor 1") +
  theme_light() + 
  annotate("text", x = 0.94, y = .12, label = purity_text)
ggsave("figures/supplementary/varComp_purity_vs_PEER_factor_1.pdf", peer1_plot, width = 4, height = 4)

#RNA concetration vs PEER 
rna_test = cor.test(purity_df$ng_ul_mean, purity_df$PEER_factor_1, use = "complete.obs")
rna_text = paste0("r=",round(rna_test$estimate, 2),", " ,"p=", signif(rna_test$p.value,3))
peer1_plot_2 = ggplot(purity_df, aes(x = ng_ul_mean, y = PEER_factor_1)) + 
  geom_point() +
  xlab("RNA concentration (ng/ul)") + 
  ylab("PEER factor 1") + 
  theme_light() + 
  annotate("text", x = 290, y = .12, label = rna_text)
ggsave("figures/supplementary/varComp_rna_concetration_vs_PEER_factor_1.pdf", peer1_plot_2, width = 4, height = 4)

#Focus on the relationship between sequencing chemistry and first PEER factor
lib_test = cor.test(as.numeric(sample_metadata_list$naive$rna_auto) , sample_metadata_list$naive$PEER_factor_2, use = "complete.obs")
lib_text = paste0("r=",round(lib_test$estimate, 2),", " ,"p<2.2e-16")
peer2_plot = ggplot(sample_metadata_list$naive, aes(x = rna_auto, y = PEER_factor_2)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = .1)) + 
  xlab("Automatic RNA-seq") + 
  ylab("PEER factor 2") +
  theme_light() +
  annotate("text", x = 1, y = .1, label = lib_text)
ggsave("figures/supplementary/varComp_autmatic_vs_PEER_factor_2.pdf", peer2_plot, width = 4, height = 4)

#### Look at the variation explained by each PEER factor ####
prec_naive = read.table("results/SL1344/PEER/naive_10/precision.txt", sep = ",")
prec_ifng = read.table("results/SL1344/PEER/IFNg_10/precision.txt", sep = ",")
prec_sl1344 = read.table("results/SL1344/PEER/SL1344_10//precision.txt", sep = ",")
prec_ifng_sl1344 = read.table("results/SL1344/PEER/IFNg_SL1344_10//precision.txt", sep = ",")

#Combine it all together
precision_df = 1/cbind(prec_naive, prec_ifng, prec_sl1344, prec_ifng_sl1344)[-1,]
colnames(precision_df) = c("naive", "IFNg", "SL1344", "IFNg_SL1344")
precision_df = dplyr::mutate(precision_df, factor = as.character(c(1:10))) %>% 
  dplyr::mutate(factor = factor(factor, factor)) %>%
  tidyr::gather(condition, sd, naive:IFNg_SL1344) %>%
  dplyr::group_by(condition) %>% 
  dplyr::mutate(var_exp = sd/sum(sd))

#Make a plot
precision_plot = ggplot(precision_df, aes(x = factor, y = var_exp, color = condition, group = condition)) +
  geom_point() + 
  geom_line() +
  ylab("Variance explained (%)") +
  theme_light() +
  xlab("PEER factor") +
  scale_color_manual(values = conditionPalette())
precision_plot
ggsave("figures/supplementary/varComp_PEER_explained.pdf", precision_plot, width = 5, height = 3.5)



