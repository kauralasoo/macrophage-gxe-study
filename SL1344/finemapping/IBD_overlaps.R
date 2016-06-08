library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
library("purrr")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#RUN MatrixEQTL on all genes
covariate_names = c("PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6")

#Extract separate lists for each condition
condition_names = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
rna_conditions = lapply(condition_names, extractConditionFromExpressionList, combined_expression_data)

#Rename column names to genotype ids
rna_conditions_renamed = lapply(rna_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#Extract gene positions
gene_pos = constructMatrixEQTLGenePos(combined_expression_data$gene_metadata)

#Extract covariates for each condition
covariates_list = purrr::map(rna_conditions_renamed, ~constructMatrixEQTLCovariates(sample_metadata = .$sample_metadata, covariate_names = covariate_names, id_column = "genotype_id"))

# RGS14
#Selected genes
selected_genes = c("ENSG00000171855","ENSG00000169220")
ifng_sl1344 = runMatrixEQTL(rna_conditions_renamed$IFNg_SL1344$cqn[selected_genes,], vcf_file$genotypes[,colnames(covariates_list$naive)], 
                  as.data.frame(vcf_file$snpspos), as.data.frame(gene_pos), covariates = as.matrix(covariates_list$IFNg_SL1344), pvOutputThreshold = 1)
ifng_sl1344_cis_qtls = matrixEQTLExtractCisQTLs(ifng_sl1344, vcf_file$snpspos)
rgs14_pvalues_ifng_sl1344 = dplyr::filter(ifng_sl1344_cis_qtls, gene_id == "ENSG00000169220")  %>%
  addR2FromLead(vcf_file$genotypes)

sl1344 = runMatrixEQTL(rna_conditions_renamed$SL1344$cqn[selected_genes,], vcf_file$genotypes[,colnames(covariates_list$naive)], 
                            as.data.frame(vcf_file$snpspos), as.data.frame(gene_pos), covariates = as.matrix(covariates_list$SL1344), pvOutputThreshold = 1)
sl1344_cis_qtls = matrixEQTLExtractCisQTLs(sl1344, vcf_file$snpspos)
rgs14_pvalues_sl1344 = dplyr::filter(sl1344_cis_qtls, gene_id == "ENSG00000169220") %>%
  addR2FromLead(vcf_file$genotypes)
ggplot(rgs14_pvalues_ifng_sl1344, aes(x = pos, y = -log(p_nominal,10), color = R2)) + geom_point()

#Look at corresponding ATAC qtls
#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Fetch corresponding SNPs from ATAC data
atac_tabix_list = list(naive = "../macrophage-chromatin/results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz",
                       IFNg = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz",
                       SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz",
                       IFNg_SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz")


rasqualTools::tabixFetchSNPsQuick("rs56235845", atac_tabix_list$IFNg, vcf_file$snpspos) %>% 
  dplyr::arrange(p_nominal) %>% dplyr::select(gene_id, p_nominal, n_feature_snps)

peak_manhattan = rasqualTools::tabixFetchGenesQuick("ATAC_peak_223859", atac_tabix_list$IFNg, atac_list$gene_metadata)[[1]] %>% 
  dplyr::arrange(p_nominal) %>%
  addR2FromLead(vcf_file$genotypes)
ggplot(peak_manhattan, aes(x = pos, y = -log(p_nominal,10), color = R2)) + geom_point()

rasqualTools::tabixFetchSNPsQuick("rs56235845", atac_tabix_list$IFNg_SL1344, vcf_file$snpspos) %>% 
  dplyr::arrange(p_nominal) %>% dplyr::select(gene_id, p_nominal, n_feature_snps)


plotEQTL("ATAC_peak_223859", "rs56235845", atac_list$cqn, vcf_file$genotypes, 
         atac_list$sample_metadata, atac_list$gene_metadata)
plotEQTL("ENSG00000169220", "rs56235845", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)

#Looks some other examples
#ITGAL - IFNg specific, but extremely weak effect
plotEQTL("ENSG00000005844", "rs11574938", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
dplyr::filter(min_pvalue_df, gene_id == "ENSG00000005844")

plotEQTL("ENSG00000108691", "rs1024611", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
dplyr::filter(min_pvalue_df, gene_id == "ENSG00000108691")

#UBQLN4 - IFNg + Salmonella specific looks unlikely
plotEQTL("ENSG00000160803", "rs4661175", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
dplyr::filter(min_pvalue_df, gene_id == "ENSG00000160803")

IFNg_SL1344_pvalues = tabixFetchGenesQuick("ENSG00000160803", "results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.sorted.txt.gz", combined_expression_data$gene_metadata, cis_window = 5e5)[[1]]  %>% 
  dplyr::mutate(condition = "IFNg_SL1344")
ggplot(IFNg_SL1344_pvalues, aes(x = pos, y = -log(p_nominal,10))) + geom_point()

#CD32 - FCGR2A - not condition speciifc and might be driven by a repeat element
dplyr::filter(min_pvalue_df, gene_id == "ENSG00000143226")
plotEQTL("ENSG00000143226", "rs4657041", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
IFNg_SL1344_pvalues = tabixFetchGenesQuick("ENSG00000143226", "results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.sorted.txt.gz", combined_expression_data$gene_metadata, cis_window = 5e5)[[1]]  %>% 
  dplyr::mutate(condition = "IFNg_SL1344")
ggplot(IFNg_SL1344_pvalues, aes(x = pos, y = -log(p_nominal,10))) + geom_point()

