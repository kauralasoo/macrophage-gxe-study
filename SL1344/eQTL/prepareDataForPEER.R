library("MatrixEQTL")
library("ggplot2")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Load expresison dataset preapred previously by processExpressionData.R script
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
gene_id_name_map = dplyr::select(expression_dataset$gene_metadata, gene_id, gene_name)

#Discard replicate samples
design = dplyr::filter(expression_dataset$design, !(donor == "fpdj")) %>% tbl_df() %>% #Remove all fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) %>% #Remove second ougl sample (dium)
  dplyr::filter(!(donor == "mijn")) #Remove mijn (wrong line from CGAP)
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#Export data for PEER
exprs_cqn_all = expression_dataset$exprs_cqn[,design$sample_id]
expressed_genes = names(which(rowMeans(exprs_cqn_all) > 0)) #Set conservative threshold to expression level
exprs_cqn = exprs_cqn_all[expressed_genes,design$sample_id]

#Export the design matrix to add sample names to the PEER data factors afterwards
saveRDS(design, "results/SL1344/PEER/input/design_matrix.rds")

#### Save only expessed genes ####
#Condition A
cond_A_design = dplyr::filter(design, condition == "A")
cond_A_exprs = t(exprs_cqn[,cond_A_design$sample_id])
write.table(cond_A_exprs, "results/SL1344/PEER/input/cond_A_exprs.txt", row.names = FALSE, col.names = FALSE, sep = ",")

#Condition B
cond_B_design = dplyr::filter(design, condition == "B")
cond_B_exprs = t(exprs_cqn[,cond_B_design$sample_id])
write.table(cond_B_exprs, "results/SL1344/PEER/input/cond_B_exprs.txt", row.names = FALSE, col.names = FALSE, sep = ",")

#Condition C
cond_C_design = dplyr::filter(design, condition == "C")
cond_C_exprs = t(exprs_cqn[,cond_C_design$sample_id])
write.table(cond_C_exprs, "results/SL1344/PEER/input/cond_C_exprs.txt", row.names = FALSE, col.names = FALSE, sep = ",")

#Condition D
cond_D_design = dplyr::filter(design, condition == "D")
cond_D_exprs = t(exprs_cqn[,cond_D_design$sample_id])
write.table(cond_D_exprs, "results/SL1344/PEER/input/cond_D_exprs.txt", row.names = FALSE, col.names = FALSE, sep = ",")

#### Save all genes ####
#Condition A
cond_A_design = dplyr::filter(design, condition == "A")
cond_A_exprs = t(exprs_cqn_all[,cond_A_design$sample_id])
write.table(cond_A_exprs, "results/SL1344/PEER/input/cond_A_exprs.all_genes.txt", row.names = FALSE, col.names = FALSE, sep = ",")

#Condition B
cond_B_design = dplyr::filter(design, condition == "B")
cond_B_exprs = t(exprs_cqn_all[,cond_B_design$sample_id])
write.table(cond_B_exprs, "results/SL1344/PEER/input/cond_B_exprs.all_genes.txt", row.names = FALSE, col.names = FALSE, sep = ",")

#Condition C
cond_C_design = dplyr::filter(design, condition == "C")
cond_C_exprs = t(exprs_cqn_all[,cond_C_design$sample_id])
write.table(cond_C_exprs, "results/SL1344/PEER/input/cond_C_exprs.all_genes.txt", row.names = FALSE, col.names = FALSE, sep = ",")

#Condition D
cond_D_design = dplyr::filter(design, condition == "D")
cond_D_exprs = t(exprs_cqn_all[,cond_D_design$sample_id])
write.table(cond_D_exprs, "results/SL1344/PEER/input/cond_D_exprs.all_genes.txt", row.names = FALSE, col.names = FALSE, sep = ",")
