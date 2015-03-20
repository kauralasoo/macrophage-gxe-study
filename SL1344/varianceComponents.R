library("dplyr")
library("lme4")
library("devtools")
library("ggplot2")
load_all("macrophage-gxe-study/seqUtils/")

varianceExplained <- function(lmer_model){
  variance = as.data.frame(VarCorr(lmer_model))
  var_percent = dplyr::mutate(variance, percent_variance = vcov/sum(vcov)) %>% 
    dplyr::select(grp, percent_variance)
  return(var_percent)  
}


#Import line metadata
line_metadata = readRDS("results/SL1344/sample_info/compiled_line_metadata.rds")
expression_list = readRDS("results/SL1344/combined_expression_data.rds")

#Remove fpdj_2 sample
design = dplyr::filter(expression_list$design, !(donor == "fpdj" & replicate == 2))
rownames(design) = design$sample_id

#Calculate mean expression per condition
exprs_cqn = expression_list$exprs_cqn[,rownames(design)]
exprs_cqn_mean = calculateMean(expression_list$exprs_cqn, expression_list$design, "condition")
expressed_genes = names(which(apply(exprs_cqn_mean, 1, max) > 1))
exprs_cqn_filtered = exprs_cqn[expressed_genes,]

#Perform PCA
pca_list = performPCA(exprs_cqn_filtered, design)
pca_plot = ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition_name, label = sample_id)) + 
  geom_text()

#Perform variance component analysis
selected_meta = dplyr::select(line_metadata, donor, replicate, max_purity, salmonella, ng_ul_mean)
gene_exp = exprs_cqn["ENSG00000183709",]
gene_df = data_frame(sample_id = names(gene_exp), exp_value = gene_exp)
model_data = dplyr::left_join(gene_df, design, by = "sample_id") %>%
  dplyr::left_join(selected_meta, by = c("donor", "replicate")) %>%
  dplyr::mutate(is_pure = ifelse(model_data$max_purity < 0.98, "no","yes"))

#Add other covariates

#Fit a model
lmer_model = lme4::lmer(exp_value ~ (1|SL1344) + (1|IFNg) + (1|IFNg:SL1344) + max_purity, model_data)
summary(lmer_model)
varianceExplained(lmer_model)

lmer_model = lme4::lmer(exp_value ~ SL1344 + IFNg + IFNg:SL1344 + ng_ul_mean + (1|donor) + (1|salmonella), model_data)

lmer_model = lme4::lmer(exp_value ~ (1|SL1344) + (1|IFNg) + (1|IFNg:SL1344) + (1|salmonella), model_data)
varianceExplained(lmer_model)

lmer_model = lme4::lmer(exp_value ~ IFNg + (IFNg|SL1344), model_data)
varianceExplained(lmer_model)

lmer_model = lme4::lmer(exp_value ~ (1|IFNg:SL1344) + (1|salmonella), model_data)
varianceExplained(lmer_model)

lm_model = lm(exp_value ~ SL1344 + IFNg + IFNg:SL1344 + factor(salmonella) + max_purity + ng_ul_mean + donor, model_data)
lm_an = anova(lm_model)
lm_an[,2] = lm_an[,2] / sum(lm_an[,2])

lmer_model = lme4::lmer(exp_value ~ (1|SL1344) + (1|IFNg) + (1|IFNg:SL1344) + (1|salmonella) + (1|donor) + ng_ul_mean, model_data)
a = varianceExplained(lmer_model)

aov_mod = aov(exp_value ~ ng_ul_mean, model_data)
a = varianceExplained(lmer_model)

library(faraway)
data(pulp)
pulp_aov = aov(bright ~ operator, pulp)

