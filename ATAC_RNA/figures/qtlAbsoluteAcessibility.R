library("devtools")
library("dplyr")
library("ggplot2")
library("purrr")
load_all("../seqUtils/")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")

#Import condition-specific QTLs
variable_qtls = readRDS("results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Only keep samples with high quality data in all conditions
hq_donors = c("qolg","vass","kuxp","cicb", "febc","eiwy","oapg","nukw","hayt","bima","pamv","guss","eipl","iill","podx","pelm")
filtered_metadata = dplyr::filter(atac_list$sample_metadata, donor %in% hq_donors)

#Calculate mean accessibility in each condition
mean_tpm = calculateMean(atac_list$tpm[,filtered_metadata$sample_id], as.data.frame(filtered_metadata), "condition_name")
mean_tpm = dplyr::mutate(mean_tpm, gene_id = rownames(mean_tpm)) %>% tbl_df()
mean_cqn = calculateMean(atac_list$cqn[,filtered_metadata$sample_id], as.data.frame(filtered_metadata), "condition_name")
mean_cqn = dplyr::mutate(mean_cqn, gene_id = rownames(mean_cqn)) %>% tbl_df()

#Extract appear QTLs
appear_qtls = variable_qtls$appear %>% 
  dplyr::select(gene_id, new_cluster_id) %>% 
  unique()

#Extract mean accessibility values
cluster_access = dplyr::filter(mean_tpm, gene_id %in% appear_qtls$gene_id) %>% 
  tidyr::gather(condition_name, tpm, IFNg:SL1344) %>% 
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
joint_access = dplyr::left_join(appear_qtls, cluster_access, by = "gene_id") %>%
  dplyr::left_join(figureNames(), by = "condition_name")

#Make a plot
caqtl_mean_accessibility = ggplot(joint_access, aes(x = figure_name, y = log(tpm + 1,2) , color = condition_name)) + 
  geom_violin(scale = "width") + 
  geom_boxplot(width = .2) +
  facet_grid(new_cluster_id~.) +
  theme_light() +
  scale_color_manual(values = conditionPalette()) +
  ylab(expression(paste("Mean accessibility (",Log[2], "(TPM+1))", sep = ""))) +
  theme(legend.position = "none", axis.title.x = element_blank())

ggsave("figures/supplementary/caQTL_appear_mean_accessibility.pdf", plot = caqtl_mean_accessibility, width = 3, height = 4.5)

##### Extract TPMs only in the naive and one stimulated condition
#Map cluster_ids to maximal conditions
max_effect_df = data_frame(new_cluster_id = c(1,2,3,4,5,6), max_effect = c("IFNg_SL1344", "IFNg_SL1344", "SL1344", "IFNg", "IFNg","IFNg"))
joint_max_df = dplyr::left_join(joint_access, max_effect_df, by = "new_cluster_id")

#Extract relevant TPMs
naive_tpms = dplyr::filter(joint_max_df, condition_name == "naive") %>% 
  dplyr::mutate(state = "Naive")
stimulated_tpms = dplyr::filter(joint_max_df, condition_name == max_effect) %>% dplyr::mutate(state = "Stimulated")
caqtl_stimulated_df = dplyr::bind_rows(naive_tpms, stimulated_tpms) %>% dplyr::mutate(phenotype = "ATAC-seq")




#### Perform the same analysis for eQTL data
#Import condition-specific QTLs
rna_variable_qtls = readRDS("results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")
rna_list = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Calculate mean accessibility in each condition
mean_tpm = calculateMean(rna_list$tpm, as.data.frame(rna_list$sample_metadata), "condition_name")
mean_tpm = dplyr::mutate(mean_tpm, gene_id = rownames(mean_tpm)) %>% tbl_df()

#Extract appear QTLs
appear_qtls = rna_variable_qtls$appear %>% 
  dplyr::ungroup() %>%
  dplyr::select(gene_id, new_cluster_id) %>% 
  unique()

#Extract mean accessibility values
cluster_access = dplyr::filter(mean_tpm, gene_id %in% appear_qtls$gene_id) %>% 
  tidyr::gather(condition_name, tpm, IFNg:SL1344) %>% 
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
joint_access = dplyr::left_join(appear_qtls, cluster_access, by = "gene_id") %>%
  dplyr::left_join(figureNames(), by = "condition_name")

##### Extract TPMs only in the naive and one stimulated condition
#Map cluster_ids to maximal conditions
max_effect_df = data_frame(new_cluster_id = c(1,2,3,4,5,6), max_effect = c("IFNg_SL1344", "IFNg_SL1344", "SL1344", "SL1344", "IFNg","IFNg"))
joint_max_df = dplyr::left_join(joint_access, max_effect_df, by = "new_cluster_id")

#Extract relevant TPMs
naive_tpms = dplyr::filter(joint_max_df, condition_name == "naive") %>% 
  dplyr::mutate(state = "Naive")
stimulated_tpms = dplyr::filter(joint_max_df, condition_name == max_effect) %>% dplyr::mutate(state = "Stimulated")
eqtl_stimulated_df = dplyr::bind_rows(naive_tpms, stimulated_tpms) %>% dplyr::mutate(phenotype = "RNA-seq")


#Make a plot
eqtl_mean_accessibility = ggplot(joint_access, aes(x = figure_name, y = log(tpm + 1,2) , color = condition_name)) + 
  geom_violin(scale = "width") + 
  geom_boxplot(width = .2) +
  facet_grid(new_cluster_id~.) +
  theme_light() +
  scale_color_manual(values = conditionPalette()) +
  ylab(expression(paste("Mean expression (",Log[2], "(TPM+1))", sep = ""))) +
  theme(legend.position = "none", axis.title.x = element_blank())
ggsave("figures/supplementary/eQTL_appear_mean_accessibility.pdf", plot = eqtl_mean_accessibility, width = 3, height = 4.5)


#Combine ATAC-seq and RNA-seq stimulated dfs together
combined_df = dplyr::bind_rows(eqtl_stimulated_df, caqtl_stimulated_df) %>%
  dplyr::mutate(phenotype = factor(phenotype, levels = c("RNA-seq", "ATAC-seq")))
combined_plot = ggplot(combined_df, aes(x = state, y = log(tpm + 1,2), color = phenotype)) + 
  geom_violin(scale = "width") + 
  geom_boxplot(width = .15, outlier.shape = NA) +
  theme_light() +
  scale_color_manual(values = c("#5e3c99","#e66101")) +
  facet_wrap(~phenotype, scale = "free_y") +
  ylab(expression(paste("Mean signal [",Log[2], "(TPM+1)]", sep = ""))) +
  xlab("Condition") +
  theme(legend.position = "none")
ggsave("figures/main_figures/mean_signal_by_phenotype.pdf", plot = combined_plot, width = 3.5, height = 3)









#Extract disappear QTLs
appear_qtls = variable_qtls$disappear %>% 
  dplyr::select(gene_id, cluster_id) %>% 
  unique()

#Extract mean accessibility values
cluster_access = dplyr::filter(mean_tpm, gene_id %in% appear_qtls$gene_id) %>% 
  tidyr::gather(condition_name, tpm, IFNg:SL1344) %>% 
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
joint_access = dplyr::left_join(appear_qtls, cluster_access, by = "gene_id")

#Make a plot
caqtl_mean_accessibility = ggplot(joint_access, aes(x = condition_name, y = log(tpm + 1,2) , color = condition_name)) + 
  geom_boxplot() + 
  facet_grid(cluster_id~.) +
  theme_light() +
  scale_color_manual(values = conditionPalette()) +
  ylab("Mean accessibility (Log2 TPM)") +
  theme(axis.text.x = element_text(angle = 15), legend.position = "none", axis.title.x = element_blank())
