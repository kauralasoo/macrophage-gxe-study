library("devtools")
load_all("macrophage-gxe-study/seqUtils/")
library(plyr)
library(dplyr)
library(ggplot2)

#Load the raw eQTL dataset
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/SL1344/array_genotypes.59_samples.vcfToMatrix.rds") #genotypes
gene_id_name_map = dplyr::select(expression_dataset$gene_metadata, gene_id, gene_name)

#Filter out some samples
design = dplyr::filter(expression_dataset$design, !(donor == "fpdj")) %>% tbl_df() %>% #Remove all fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) #Remove second ougl sample (dium)
selected_samples = design$sample_id
expression_dataset$design = expression_dataset$design[selected_samples,]
expression_dataset$exprs_cqn = expression_dataset$exprs_cqn[,selected_samples]
expression_dataset$exprs_counts = expression_dataset$exprs_counts[,selected_samples]

#Load identified eQTLS
condA_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_A.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map, fdr_cutoff = 0.1)
condB_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_B.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map, fdr_cutoff = 0.1)
condC_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_C.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map, fdr_cutoff = 0.1)
condD_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_D.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map, fdr_cutoff = 0.1)

#For each gene keep the eqtl with lowest p-value
joint_eqtls = rbind(condA_eqtls, condB_eqtls, condC_eqtls, condD_eqtls) %>%
  dplyr::group_by(gene_id) %>% 
  arrange(pvalue) %>% 
  dplyr::filter(row_number() == 1)

#Test interaction between genotype and condition
res = testMultipleInteractions(joint_eqtls, expression_dataset, vcf_file, line_metadata)
gene_id_name_map = dplyr::select(expression_dataset$gene_metadata, gene_id, gene_name
interaction_pvalues = ldply(res,function(x){x[[6]][2]},.id = "gene_id") %>%
  dplyr::arrange(V1)
interaction_df = dplyr::left_join(interaction_pvalues, joint_eqtls, by = "gene_id")
interaction_df = dplyr::filter(interaction_df,V1 < 0.05)

interaction_plots = makeMultiplePlots(interaction_df, expression_dataset, vcf_file, line_metadata)
savePlots(interaction_plots, path = "results/SL1344/matrixeQTL/interaction_plots/",width = 7, height = 7)

#eQTLS found in SL1344
condA_plots = makeMultiplePlots(condA_eqtls, expression_dataset, vcf_file, line_metadata)
savePlots(condA_plots, path = "results/SL1344/matrixeQTL/condA_plots/",width = 7, height = 7)

#eQTLS found in IFNg
condB_plots = makeMultiplePlots(condB_eqtls, expression_dataset, vcf_file, line_metadata)
savePlots(condB_plots, path = "results/SL1344/matrixeQTL/condB_plots/",width = 7, height = 7)

#eQTLS found in SL1344
condC_plots = makeMultiplePlots(condC_eqtls, expression_dataset, vcf_file, line_metadata)
savePlots(condC_plots, path = "results/SL1344/matrixeQTL/condC_plots/",width = 7, height = 7)

#eQTLS found in SL1344
condD_plots = makeMultiplePlots(condD_eqtls, expression_dataset, vcf_file, line_metadata)
savePlots(condD_plots, path = "results/SL1344/matrixeQTL/condD_plots/",width = 7, height = 7)


condA_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_A.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
condB_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_B.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
condC_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_C.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
condD_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_D.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
eqtl_list = list(A = condA_eqtls, B = condB_eqtls, C = condC_eqtls, D = condD_eqtls)




