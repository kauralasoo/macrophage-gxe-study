library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")
library("readr")
load_all("~/software/rasqual/rasqualTools/")

#Import atac data list
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Extract separate lists for each condition
naive_list = extractConditionFromExpressionList("naive", atac_list)
IFNg_list = extractConditionFromExpressionList("IFNg", atac_list)
SL1344_list = extractConditionFromExpressionList("SL1344", atac_list)
IFNg_SL1344_list = extractConditionFromExpressionList("IFNg_SL1344", atac_list)
atac_conditions = list(naive = naive_list, IFNg = IFNg_list, SL1344 = SL1344_list, IFNg_SL1344 = IFNg_SL1344_list)

#Save expression data for PEER and run PEER outside of R
peer_cqn_list = lapply(atac_conditions, function(x){x$cqn})
savePEERData(peer_cqn_list, "results/ATAC/PEER/input/")

peer_tpm_list = lapply(atac_conditions, function(x){log(x$tpm + 0.01, 2)})
savePEERData(peer_cqn_list, "results/ATAC/PEER/input/",file_suffix = "exprs_tpm")

#Import PEER factors back into R
naive_peer = importPEERFactors("results/ATAC/PEER/output/naive_10/factors.txt", atac_conditions$naive$sample_metadata)
IFNg_peer = importPEERFactors("results/ATAC/PEER/output/IFNg_10/factors.txt", atac_conditions$IFNg$sample_metadata)
SL1344_peer = importPEERFactors("results/ATAC/PEER/output/SL1344_10/factors.txt", atac_conditions$SL1344$sample_metadata)
IFNg_SL1344_peer = importPEERFactors("results/ATAC/PEER/output/IFNg_SL1344_10/factors.txt", atac_conditions$IFNg_SL1344$sample_metadata)
peer_factors = rbind(naive_peer, IFNg_peer, SL1344_peer, IFNg_SL1344_peer)
sample_ids = peer_factors$sample_id
covariates = t(peer_factors[,-1])
colnames(covariates)= sample_ids
covariates = covariates[,colnames(head(atac_list$counts))]

#Use PCA to construct covariates for QTL mapping
naive_PCA = performPCA(log(atac_conditions$naive$tpm + 0.01, 2), atac_conditions$naive$sample_metadata, n_pcs = 10)
IFNg_PCA = performPCA(log(atac_conditions$IFNg$tpm + 0.01, 2), atac_conditions$IFNg$sample_metadata, n_pcs = 10)
SL1344_PCA = performPCA(log(atac_conditions$SL1344$tpm + 0.01, 2), atac_conditions$SL1344$sample_metadata, n_pcs = 10)
IFNg_SL1344_PCA = performPCA(log(atac_conditions$IFNg_SL1344$tpm + 0.01, 2), atac_conditions$IFNg_SL1344$sample_metadata, n_pcs = 10)

#Plot variance explained
plot(naive_PCA$var_exp)
plot(IFNg_PCA$var_exp)
plot(SL1344_PCA$var_exp)
plot(IFNg_SL1344_PCA$var_exp)

#Add PCs into the sample metadata df
covariates = rbind(naive_PCA$pca_matrix, IFNg_PCA$pca_matrix, SL1344_PCA$pca_matrix, IFNg_SL1344_PCA$pca_matrix) %>% 
  tbl_df() %>%
  dplyr::mutate(sex_binary = ifelse(sex == "male",0,1))
new_sample_metadata = dplyr::select(atac_list$sample_metadata, sample_id) %>% dplyr::left_join(covariates, by = "sample_id")
atac_list$sample_metadata = new_sample_metadata
