library("plyr")
library("dplyr")
library("lme4")
library("devtools")
library("ggplot2")
load_all("macrophage-gxe-study/seqUtils/")

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

#Specify alternative models
model1 <- function(model_data){
  mod = lme4::lmer(exp_value ~ (1|SL1344) + (1|IFNg) + (1|SL1344:IFNg) + (1|salmonella) + (1|donor), model_data)
  return(mod)
}

#Perform variance component analysis
selected_meta = dplyr::select(line_metadata, donor, replicate, max_purity, salmonella, ng_ul_mean)

#Apply lmer to a list of genes
gene_ids = rownames(exprs_cqn_filtered)
gene_id_list = idVectorToList(gene_ids)
gene_data_list = lapply(gene_id_list, constructGeneData, exprs_cqn, design, selected_meta)
variance_list = lapply(gene_data_list, estimateVarianceExplained, model1)
var_table = ldply(variance_list, .id = "gene_id")
saveRDS(var_table, "results/varComp/model1_results.rds")

#Plot binned estimates of variance explained
var_table = readRDS("results/varComp/model1_results.rds") %>% tbl_df()

binned_table = dplyr::arrange(var_table, Residual) %>% 
  dplyr::mutate(residual_bin = 20-floor(Residual*20)) %>% 
  dplyr::select(-type) %>% 
  dplyr::rename(SL1344_IFNg = `SL1344:IFNg`) %>%
  dplyr::rename(batch = salmonella) %>%
  dplyr::rename(line = donor)

#Calculate the number of genes in each bin
bin_size_table = dplyr::group_by(binned_table, residual_bin) %>% dplyr::summarise(bin_size = length(residual_bin))

#Calculate mean variance explained within each group
var_explained = tidyr::gather(binned_table, component, var_explained, line:SL1344_IFNg) %>% 
  group_by(residual_bin, component) %>% 
  dplyr::summarise(var_explained = mean(var_explained)) %>%
  dplyr::mutate(component = factor(as.vector(component), 
                  levels = c("Residual", "SL1344", "IFNg", "SL1344_IFNg", "line", "batch"))) %>%
  dplyr::arrange(component)

#Make plots
bin_sizes_plot = ggplot(bin_size_table, aes(x = residual_bin, y = bin_size)) + geom_bar(stat = "identity") + 
  ylab("Number of genes") + 
  xlab("Residual variance bin")
ggsave("results/varComp/bin_sizes_plot.pdf", plot = bin_sizes_plot, width = 8, height = 3)
var_comp_plot = ggplot(var_explained, aes(x = residual_bin, y = var_explained, fill = component)) + geom_bar(stat="identity") +
  ylab("% variance explained") +
  xlab("Residual variance bin")
ggsave("results/varComp/var_comp_plot.pdf", plot = var_comp_plot, width = 11, height = 7)


#Fit a model
summary(lmer_model)
varianceExplained(lmer_model)

summary(lmer_model)
varianceExplained(lmer_model)


lmer_model = lme4::lmer(exp_value ~ SL1344 + IFNg + IFNg:SL1344 + (1|donor) + (1|salmonella), model_data)

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

warningTest <- function(x){
  if(x < 5){
    warning("X is too small.")
  }
  return(x)
}

run <- function(x){
  tryCatch({
    d = warningTest(x)
    return(list(value = d, converged = TRUE))
    },
    warning = function(c){
      d = suppressWarnings(warningTest(x))
      return(list(value = d, converged = FALSE))
    })
}
