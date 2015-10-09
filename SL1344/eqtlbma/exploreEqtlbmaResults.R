library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("ggplot2")
load_all("macrophage-gxe-study/housekeeping/")
library("tidyr")
library("lme4")

#Load the raw eQTL dataset
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/SL1344/array_genotypes.59_samples.imputed.vcfToMatrix.rds") #genotypes
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Filter the expression data set
selected_samples = eqtl_data_list$sample_metadata$sample_id
expression_dataset$design = expression_dataset$design[selected_samples,]
expression_dataset$exprs_cqn = expression_dataset$exprs_cqn[,selected_samples]
expression_dataset$exprs_counts = expression_dataset$exprs_counts[,selected_samples]

#Import eqtlbma results
results = read.table(gzfile("eqtlbma/output_imputed//eqtlbma_imputed_avg_bfs.txt.gz"), sep = "\t", skip = 1, header = TRUE, stringsAsFactors = FALSE)
significant_genes = dplyr::filter(results, gene.post > 0.8) %>% tbl_df() %>%
  dplyr::select(gene, snp, gene.post, snp.post.the, snp.post.an, best.config, post.best.config) %>%
  dplyr::group_by(gene) %>%
  dplyr::arrange(-snp.post.an) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::rename(gene_id = gene, snp_id = snp) %>%
  dplyr::left_join(gene_id_name_map, by = "gene_id")
table(significant_genes$best.config)

#Make plots for condition-specific eQTLs
dplyr::filter(significant_genes, best.config == "1-2-3") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
savePlots(.,path = "results/SL1344/eqtlbma/inflammatory_plots/",width = 7, height = 7)

#Make plots for condition-specific eQTLs
dplyr::filter(significant_genes, best.config %in% c("1-2","1")) %>%
  makeMultiplePlots(., expression_dataset, vcf_file, line_metadata) %>%
  savePlots(., path = "results/SL1344/eqtlbma/ifng_plots/",width = 7, height = 7)

#Salmonella specific
dplyr::filter(significant_genes, best.config == "2-3") %>%
  makeMultiplePlots(., expression_dataset, vcf_file, line_metadata) %>%
  savePlots(., path = "results/SL1344/eqtlbma/SL1344_plots/",width = 7, height = 7)

# Not Salmonella specific
dplyr::filter(significant_genes, best.config %in% c("1-4","1-3-4")) %>%
  makeMultiplePlots(., expression_dataset, vcf_file, line_metadata) %>%
  savePlots(., path = "results/SL1344/eqtlbma/not_SL1344_plots/",width = 7, height = 7)

#Test for interactions between lead SNP and condition
res = testMultipleInteractions(significant_genes, eqtl_data_list)
interaction_pvalues = ldply(res,function(x){x$anova[[6]][2]},.id = "gene_id") %>%
  dplyr::arrange(V1) %>%
  dplyr::rename(interaction = V1)
interaction_df = dplyr::left_join(interaction_pvalues, significant_genes, by = "gene_id") %>%
  dplyr::mutate(interaction_fdr = p.adjust(interaction, "fdr"))

interaction_qtls = dplyr::filter(interaction_df,interaction_fdr < 0.1)
interaction_plots = makeMultiplePlots(interaction_qtls, expression_dataset, vcf_file, line_metadata)
savePlots(interaction_plots, path = "results/SL1344/eqtlbma/interaction_plots/",width = 7, height = 7)
write.table(interaction_qtls, "results/SL1344/interaction_qtl_results.txt", quote = FALSE, sep = "\t")

#Explore the raw BFs and not posteriors
a = results %>% tbl_df() %>% 
  dplyr::rename(gene_id = gene, snp_id = snp) %>% 
  dplyr::semi_join(significant_genes, by = c("gene_id"))
filtered_results = dplyr::semi_join(a, significant_genes, by = c("gene_id","snp_id")) %>% 
  dplyr::select(gene_id, snp_id, log10.bf.1.2.3, log10.bf.1.2, log10.bf.2.3,log10.bf.1.2.3.4,log10.bf.1.4,log10.bf.3.4, log10.bf.4,log10.bf.2,log10.bf.1.3.4) %>%
  dplyr::left_join(gene_id_name_map, by = "gene_id") %>%
  dplyr::select(gene_id, gene_name, snp_id, everything()) %>%
  dplyr::semi_join(interaction_qtls, by = "gene_id")

#Find the best configuration for each gene
config_labels = data_frame(best_config_index = c(1,2,3,4,5,6,7,8,9), 
                           best_config = c("Both", "IFNg", "Salmonella", "Shared", "not_Salmonella", "not_IFNg","not_Both","Synergistic","not_Synergistic"))
best_config = dplyr::mutate(filtered_results, best_config_index = apply(filtered_results[,4:12],1, which.max)) %>%
  dplyr::left_join(config_labels, by = "best_config_index")

#Count different types of eQTLs
config_counts = dplyr::group_by(best_config, best_config) %>% 
  dplyr::summarise(config_count = length(best_config)) %>% 
  dplyr::mutate(mode = c(rep("appeared",5), rep("disappeared",4))) %>% 
  dplyr::mutate(mode = ifelse(best_config == "Shared", "shared", mode))
config_count_plot = ggplot(config_counts, aes(x = best_config, y = config_count)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~mode, scales = "free_x", space = "free_x")
ggsave("results/SL1344/eqtlbma/config_counts.png",plot = config_count_plot, width = 8, height = 6)

#Make plots for each config
dplyr::filter(best_config, best_config == "Inflammatory") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
  savePlots(.,path = "results/SL1344/eqtlbma/two-step/Inflammatory/",width = 7, height = 7)

dplyr::filter(best_config, best_config == "IFNg") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
  savePlots(.,path = "results/SL1344/eqtlbma/two-step/IFNg/",width = 7, height = 7)

dplyr::filter(best_config, best_config == "SL1344") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
  savePlots(.,path = "results/SL1344/eqtlbma/two-step/SL1344/",width = 7, height = 7)

dplyr::filter(best_config, best_config == "Shared") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
  savePlots(.,path = "results/SL1344/eqtlbma/two-step/Shared/",width = 7, height = 7)

dplyr::filter(best_config, best_config == "not_SL1344") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
  savePlots(.,path = "results/SL1344/eqtlbma/two-step/not_SL1344/",width = 7, height = 7)

dplyr::filter(best_config, best_config == "not_IFNg") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
  savePlots(.,path = "results/SL1344/eqtlbma/two-step/not_IFNg/",width = 7, height = 7)

dplyr::filter(best_config, best_config == "naive") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
  savePlots(.,path = "results/SL1344/eqtlbma/two-step/naive/",width = 7, height = 7)

dplyr::filter(best_config, best_config == "IFNg_SL1344") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
  savePlots(.,path = "results/SL1344/eqtlbma/two-step/IFNg_SL1344/",width = 7, height = 7)

dplyr::filter(best_config, best_config == "not_IFNg_SL1344") %>%
  makeMultiplePlots(.,expression_dataset, vcf_file, line_metadata) %>%
  savePlots(.,path = "results/SL1344/eqtlbma/two-step/not_IFNg_SL1344/",width = 7, height = 7)


#Play around with lmer4
golga7_data = constructGeneData("ENSG00000147533","rs6998193", eqtl_data_list)
lmer2_model = lmer(expression ~  genotype + (1+genotype|condition_name), golga7_data, REML = FALSE)
lmer1_model = lmer(expression ~  genotype + (1|condition_name), golga7_data, REML = FALSE)
a = anova(lmer1_model, lmer2_model)


golga7_data = constructGeneData("ENSG00000170458","rs778583", eqtl_data_list)
lmer2_model = lmer(expression ~  genotype + (1+genotype|condition_name), golga7_data, REML = FALSE)
lmer1_model = lmer(expression ~  genotype + (1|condition_name), golga7_data, REML = FALSE)
a = anova(lmer1_model, lmer2_model)

lmer3_model = lmer(expression ~ (1|condition_name), golga7_data, REML = FALSE)

lmer1_model = lmer(expression ~  genotype + (1+genotype|condition_name), golga7_data, REML = FALSE)
lmer2_model = lmer(expression ~  (1+genotype|condition_name), golga7_data, REML = FALSE)
anova(lmer1_model, lmer2_model)

lm1_model = lm(expression ~  genotype, golga7_data)
lm2_model = lm(expression ~ 1, golga7_data)
anova(lm1_model, lm2_model)