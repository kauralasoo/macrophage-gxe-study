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
dat = tidyr::gather(var_table, factor, var_explained, chemistry:IFNg_SL1344)
variance_dist_plot = ggplot(dat, aes(x = factor, y = var_explained)) + 
  geom_violin(scale = "width" ) + 
  geom_boxplot(width = .2) + 
  theme_light()
variance_dist_plot
ggsave("figures/supplementary/varComp_variance_dist_plot.pdf", variance_dist_plot, width = 10, height = 7)




