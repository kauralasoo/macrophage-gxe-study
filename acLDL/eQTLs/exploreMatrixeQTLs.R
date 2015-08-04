library("ggplot2")
library("devtools")
library("plyr")
library("dplyr")
load_all("macrophage-gxe-study/seqUtils/")

#Load expression dataset
dataset = readRDS("results/acLDL/acLDL_combined_expression_data.rds")
dataset$design$replicate = 1 #Add replicate to the desig
colnames(dataset$design)[4] = "condition_name"
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/acLDL/acLDL_array_genotypes.GRCh38.vcfToMatrix.rds") #genotypes
gene_id_name_map = dplyr::select(dataset$gene_metadata, gene_id, gene_name)

#Load identified eQTLS
ctrl_eqtls = filterEQTLs(read.table("results/acLDL/matrixeQTL/matrixeQTL_Ctrl.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map)
acldl_eqtls = filterEQTLs(read.table("results/acLDL/matrixeQTL/matrixeQTL_AcLDL.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map)

#Put all eQTLs together
joint_eqtls = rbind(ctrl_eqtls, acldl_eqtls) %>%
  dplyr::group_by(gene_id) %>% 
  arrange(pvalue) %>% 
  dplyr::filter(row_number() == 1)

#Remove outlier samples xegx and nusw
design = dplyr::filter(dataset$design, !(donor %in% c("xegx")))
selected_samples = design$sample_id
expressed_genes = names(which(rowMeans(dataset$exprs_cqn) > 0)) #Set conservative threshold to expression level
expressed_genes = union(joint_eqtls$gene_id, expressed_genes)

#Filter expression dataset to keep only good samples and expressed genes
dataset_filtered = filterExpressionDataset(dataset, sample_ids = selected_samples, gene_ids = expressed_genes)

#Regress outs PCs condition-by-condition
ctrl_data = dataset_filtered$exprs_cqn[,dplyr::filter(design, condition_name == "Ctrl")$sample_id]
acldl_data = dataset_filtered$exprs_cqn[,dplyr::filter(design, condition_name == "AcLDL")$sample_id]
ctrl_data_pca = regressPrinicpalComponents(ctrl_data, 4)
acldl_data_pca = regressPrinicpalComponents(acldl_data, 4)
cqn_pca  =cbind(acldl_data_pca, ctrl_data_pca)[,colnames(dataset_filtered$exprs_cqn)]
dataset_filtered$exprs_cqn = cqn_pca

#Test interaction between genotype and condition
res = testMultipleInteractions(joint_eqtls, dataset_filtered, vcf_file, line_metadata)
gene_id_name_map = dplyr::select(dataset_filtered$gene_metadata, gene_id, gene_name)
interaction_pvalues = ldply(res,function(x){x[[6]][2]},.id = "gene_id") %>%
  dplyr::arrange(V1) %>%
  dplyr::rename(interaction_pvalue = V1)
interaction_df = dplyr::left_join(interaction_pvalues, joint_eqtls, by = "gene_id") %>%
  dplyr::mutate(interaction_FDR = p.adjust(interaction_pvalue, method = "fdr"))
interaction_df = dplyr::filter(interaction_df,interaction_pvalue < 0.05)
write.table(interaction_df, "results/acLDL/matrixeQTL/interaction_eqtls.txt", sep = "\t", quote =FALSE, row.names = FALSE)

#Make plots of interaction genes
interaction_plots = makeMultiplePlots_acLDL(interaction_df, dataset_filtered, vcf_file, line_metadata)
savePlots(interaction_plots, path = "results/acLDL/matrixeQTL/interaction_plots/",width = 7, height = 7)


#eQTLS found in SL1344
ctrl_plots = makeMultiplePlots_acLDL(ctrl_eqtls, dataset_f, vcf_file, line_metadata)
savePlots(ctrl_plots, path = "results/acLDL/matrixeQTL/ctrl_plots/",width = 7, height = 7)

#eQTLS found in SL1344
acldl_plots = makeMultiplePlots_acLDL(acldl_eqtls, dataset, vcf_file, line_metadata)
savePlots(acldl_plots, path = "results/acLDL/matrixeQTL/acldl_plots/",width = 7, height = 7)

