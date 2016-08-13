library("readr")

#Import variance component analysis results
colnames = colnames(var_table)
colnames = c("gene_id","type","chemistry","RNA_concentration","IFNg","library_pool","line_id",
             "diff_duration","passage","purity",
             "residual","library_type","RNA_ext_date","stimulation_date",
             "sex","SL1344","IFNg_SL1344","converged") 
var_table = readr::read_delim("results/SL1344/varComp/varComp_all.txt", delim = "\t", col_names = colnames) %>%
  dplyr::filter(converged == TRUE) %>% #Remove genes where lmr4 did not converge
  dplyr::select(-type, -converged) #Remove unnecessary columns

#Bin genes by residual variance
binned_table = binGenesByResidual(var_table, n_bins = 20)
bin_size_table = dplyr::group_by(binned_table, residual_bin) %>% 
  dplyr::summarise(bin_size = length(residual_bin))
bin_sizes_plot = ggplot(bin_size_table, aes(x = residual_bin, y = bin_size)) + geom_bar(stat = "identity") + 
  ylab("Number of genes") + 
  xlab("Residual variance bin") + 
  theme_light()
bin_sizes_plot
ggsave("figures/supplementary/varComp_variance_bins_hist.pdf", bin_sizes_plot, width = 6, height = 2)

#Look what explains most variance within each bin
var_explained = meanVarianceWithinBins(binned_table)

#Lets look at total variance explained by each factor
total_var = dplyr::filter(var_explained, residual_bin == "Total") %>% 
  arrange(-var_explained)
total_var
write.table(total_var, "figures/supplementary/varComp_total_variance_explained.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#Keep only factors that explain more than 2% of the total variation
selected_components = as.vector(dplyr::filter(total_var, var_explained > 0.02)$component)
var_exp_selected = dplyr::filter(var_explained, component %in% selected_components) %>%
  dplyr::mutate(component = as.character(component))
plotBinnedVariance(var_exp_selected)

#We can look at the distribution of variance explained by each factor
dat = tidyr::gather(var_table, factor, var_explained, chemistry:IFNg_SL1344) %>%
  dplyr::mutate(factor = factor(factor, levels = as.character(total_var$component)))
variance_dist_plot = ggplot(dat, aes(x = factor, y = var_explained)) + 
  geom_violin(scale = "width" ) + 
  geom_boxplot(width = .2) + 
  theme_light() + 
  theme(axis.text.x=element_text(angle=15))
variance_dist_plot
ggsave("figures/supplementary/varComp_variance_dist_plot.pdf", variance_dist_plot, width = 10, height = 7)


#### Variannce components by condition ####
colnames = c("condition_name","gene_id","type","chemistry","RNA_concentration","library_pool",
             "diff_duration","passage","purity","residual","library_type",
             "RNA_ext_date","stimulation_date","sex","converged")
var_table_conditions = readr::read_delim("results/SL1344/varComp/varCompByCondition_all.txt", 
                              delim = "\t", col_names = colnames) %>%
  dplyr::filter(converged == TRUE) %>% #Remove genes where lmr4 did not converge
  dplyr::select(-type, -converged) #Remove unnecessary columns

#Gather different components together
var_table_tidy = tidyr::gather(var_table_conditions, factor, var_explained, chemistry:sex)

#Look at mean variance explained by eah factor
mean_var = dplyr::group_by(var_table_tidy, condition_name, factor) %>%
  dplyr::summarise(var_mean = mean(var_explained), var_sd = sd(var_explained)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(condition_name, -var_mean) %>%
  dplyr::mutate(ymin = var_mean - var_sd, ymax = var_mean + var_sd) %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = rev(c("naive","IFNg","SL1344", "IFNg_SL1344"))))

#Arrange factors by the amount of variance explained
factor_order = rev(as.character(dplyr::filter(mean_var, condition_name == "naive")$factor))
mean_var = dplyr::mutate(mean_var, factor = factor(factor, levels = factor_order))

#Make a line-plot of mean variance explained by each factor
mean_var_plot = dplyr::filter(mean_var, factor != "residual") %>%
  ggplot(aes(x = factor, y = var_mean, color = condition_name, group = condition_name,
             ymin = ymin, ymax = ymax)) + 
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) + 
  coord_flip() + 
  scale_color_manual(values = rev(conditionPalette())) + 
  theme_light() +
  theme(legend.key = element_blank(), axis.title.y = element_blank()) + 
  ylab("Variance explained")
ggsave("figures/supplementary/varComp_mean_var_exp_by_condition.pdf",plot = mean_var_plot, height = 6, width = 6)

#Make violin plots of the variance explained by stimulation date
stimulation_date_plot = dplyr::filter(var_table_tidy, factor == "stimulation_date") %>% 
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344", "IFNg_SL1344"))) %>%
  ggplot(aes(x = condition_name, y = var_explained, color = condition_name)) + 
  geom_violin(scale = "width" ) + 
  geom_boxplot(width = .1) + 
  theme_light() +
  scale_color_manual(values = conditionPalette()) +
  theme(legend.key = element_blank()) + 
  ylab("Variance explained") + 
  xlab("Condition")
ggsave("figures/supplementary/varComp_stimulation_date.pdf",plot = stimulation_date_plot, height = 6, width = 6)

  
  

  

