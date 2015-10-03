library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("ggplot2")

#Load the raw eQTL dataset
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/SL1344/array_genotypes.59_samples.vcfToMatrix.rds") #genotypes
gene_id_name_map = dplyr::select(expression_dataset$gene_metadata, gene_id, gene_name)

#Filter out some samples
design = dplyr::filter(expression_dataset$design, !(donor == "fpdj")) %>% tbl_df() %>% #Remove all fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) %>% #Remove second ougl sample (dium)
  dplyr::filter(!(donor == "mijn")) #%>% #Remove mijn (wrong line from CGAP)

selected_samples = design$sample_id
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

#Salmonella specific
not_SL1344_qtls = dplyr::filter(significant_genes, best.config == "1-4")
not_SL1344_plots = makeMultiplePlots(not_SL1344_qtls, expression_dataset, vcf_file, line_metadata)
savePlots(not_SL1344_plots, path = "results/SL1344/eqtlbma/not_SL1344_plots/",width = 7, height = 7)

shared_top50 = dplyr::filter(significant_genes, best.config == "1-2-3-4") %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(-gene.post) %>% 
  head(50)
shared_plots = makeMultiplePlots(shared_top50, expression_dataset, vcf_file, line_metadata)
savePlots(shared_plots, path = "results/SL1344/eqtlbma/shared_plots/",width = 7, height = 7)

