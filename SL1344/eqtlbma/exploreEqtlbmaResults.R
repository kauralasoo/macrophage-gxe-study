library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("ggplot2")
load_all("macrophage-gxe-study/housekeeping/")


#Load the raw eQTL dataset
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/SL1344/array_genotypes.59_samples.vcfToMatrix.rds") #genotypes
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Filter the expression data set
selected_samples = eqtl_data_list$sample_metadata$sample_id
expression_dataset$design = expression_dataset$design[selected_samples,]
expression_dataset$exprs_cqn = expression_dataset$exprs_cqn[,selected_samples]
expression_dataset$exprs_counts = expression_dataset$exprs_counts[,selected_samples]

#Import eqtlbma results
results = read.table(gzfile("eqtlbma/output/eqtlbma_batched_avg_bfs.txt.gz"), sep = "\t", skip = 1, header = TRUE, stringsAsFactors = FALSE)
significant_genes = dplyr::filter(results, gene.post > 0.8) %>% tbl_df() %>%
  dplyr::select(gene, snp, gene.post, snp.post.the, snp.post.an, best.config, post.best.config) %>%
  dplyr::group_by(gene) %>%
  dplyr::arrange(-snp.post.an) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::rename(gene_id = gene, snp_id = snp) %>%
  dplyr::left_join(gene_id_name_map, by = "gene_id")
table(significant_genes$best.config)

#Make plots for condition-specific eQTLs
inflammatory_qtls = dplyr::filter(significant_genes, best.config == "1-2-3")
inflammatory_plots = makeMultiplePlots(inflammatory_qtls, expression_dataset, vcf_file, line_metadata)
savePlots(inflammatory_plots, path = "results/SL1344/eqtlbma/inflammatory_plots/",width = 7, height = 7)

#Make plots for condition-specific eQTLs
ifng_qtls = dplyr::filter(significant_genes, best.config %in% c("1-2","1"))
ifng_plots = makeMultiplePlots(ifng_qtls, expression_dataset, vcf_file, line_metadata)
savePlots(ifng_plots, path = "results/SL1344/eqtlbma/ifng_plots/",width = 7, height = 7)

#Salmonella specific
SL1344_qtls = dplyr::filter(significant_genes, best.config == "2-3")
SL1344_plots = makeMultiplePlots(SL1344_qtls, expression_dataset, vcf_file, line_metadata)
savePlots(SL1344_plots, path = "results/SL1344/eqtlbma/SL1344_plots/",width = 7, height = 7)

# Not Salmonella specific
not_SL1344_qtls = dplyr::filter(significant_genes, best.config %in% c("1-4","1-3-4"))
not_SL1344_plots = makeMultiplePlots(not_SL1344_qtls, expression_dataset, vcf_file, line_metadata)
savePlots(not_SL1344_plots, path = "results/SL1344/eqtlbma/not_SL1344_plots/",width = 7, height = 7)

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


#Extract genotype effect sizes
interaction_coeficients = ldply(res,function(x){x$interaction_model$coefficients[c(2,13,14,15)]},.id = "gene_id")
rownames(interaction_coeficients) = interaction_coeficients$gene_id
interaction_coeficients = interaction_coeficients[,-1]
interaction_coeficients = cbind(interaction_coeficients[,1],interaction_coeficients[,2:4]+interaction_coeficients[,1])


interaction_qtls = dplyr::filter(interaction_df,interaction_fdr < 0.001)
selected = interaction_coeficients[interaction_qtls$gene_id,]
rownames(selected) = interaction_qtls$gene_name
heatmap.2(cor(t(zScoreNormalize(abs(selected))), method = "pearson"), tracecol = NA)

