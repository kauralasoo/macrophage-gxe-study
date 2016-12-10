library("devtools")
library("dplyr")
load_all("../seqUtils/")
load_all("../macrophage-gxe-study/macrophage-gxe-study/housekeeping/")
library("ggplot2")
library("purrr")

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
joint_access = dplyr::left_join(appear_qtls, cluster_access, by = "gene_id")

#Make a plot
caqtl_mean_accessibility = ggplot(joint_access, aes(x = condition_name, y = log(tpm + 1,2) , color = condition_name)) + 
  geom_violin(scale = "width") + 
  geom_boxplot(width = .2) +
  facet_grid(new_cluster_id~.) +
  theme_light() +
  scale_color_manual(values = conditionPalette()) +
  ylab(expression(paste("Mean accessibility (",Log[2], "(TPM+1))", sep = ""))) +
  theme(axis.text.x = element_text(angle = 15), legend.position = "none", axis.title.x = element_blank())

ggsave("figures/supplementary/caQTL_appear_mean_accessibility.pdf", plot = caqtl_mean_accessibility, width = 4, height = 7)


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
joint_access = dplyr::left_join(appear_qtls, cluster_access, by = "gene_id")

#Make a plot
eqtl_mean_accessibility = ggplot(joint_access, aes(x = condition_name, y = log(tpm + 1,2) , color = condition_name)) + 
  geom_violin(scale = "width") + 
  geom_boxplot(width = .2) +
  facet_grid(new_cluster_id~.) +
  theme_light() +
  scale_color_manual(values = conditionPalette()) +
  ylab(expression(paste("Mean expression (",Log[2], "(TPM+1))", sep = ""))) +
  theme(axis.text.x = element_text(angle = 15), legend.position = "none", axis.title.x = element_blank())

ggsave("figures/supplementary/eQTL_appear_mean_accessibility.pdf", plot = eqtl_mean_accessibility, width = 4.3, height = 7)








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
