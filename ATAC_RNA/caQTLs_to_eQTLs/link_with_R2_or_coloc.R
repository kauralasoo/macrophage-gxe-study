library("devtools")
library("plyr")
library("dplyr")
library("purrr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

#Helper functions
fetchRasqualSNPs <- function(snp_ids, snpspos, summary_list){
  #Construct GRanges object of SNP positions
  selected_snps = dplyr::filter(snpspos, snpid %in% snp_ids) %>%
    dplyr::transmute(snp_id = snpid, seqnames = chr, start = pos, end = pos, strand = "*") %>%
    dataFrameToGRanges()

  selected_pvalues = lapply(summary_list, function(tabix, snps) rasqualTools::tabixFetchSNPs(snps, tabix), selected_snps)
  return(selected_pvalues)
}

#Exptract cluster-specific gene-peak pairs from all pairs
extractGenePeakPairs <- function(all_pairs, cluter_genes){
  selected_pairs = dplyr::semi_join(all_pairs, cluter_genes, by = c("gene_id", "snp_id")) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::arrange(gene_id, -R2, p_nominal) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene_id, snp_id, peak_id) %>% 
    unique()
  return(selected_pairs)
}

#Exptract cluster-specific peak-gene pairs from all pairs (reverse analysis)
extractPeakGenePairs <- function(all_pairs, cluter_genes){
  selected_pairs = dplyr::rename(cluter_genes, peak_id = gene_id) %>%
    dplyr::semi_join(all_pairs, ., by = c("peak_id", "snp_id")) %>%
    dplyr::group_by(peak_id) %>%
    dplyr::arrange(peak_id, -R2, p_nominal) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene_id, snp_id, peak_id) %>% 
    unique()
  return(selected_pairs)
}

quantileNormaliseBeta <- function(beta){
  m = mean(beta)
  new_beta = quantileNormaliseVector(beta - m) + m
}

sortByBeta <- function(beta_df, phenotype_name){
  sorted_names = dplyr::filter(beta_df, phenotype == phenotype_name) %>% 
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(beta_mean = mean(beta)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(condition_name == "naive") %>%
    dplyr::arrange(beta_quantile)
  sorted_df = dplyr::mutate(beta_df, gene_name = factor(as.character(gene_name), 
                                                        levels = as.character(sorted_names$gene_name)))
  return(sorted_df)
}

plotQTLBetas <- function(beta_df){
  n_pairs = nrow(dplyr::select(beta_df, gene_name, snp_id) %>% unique())
  ylabel = paste(n_pairs, "eQTL-caQTL pairs")
  effect_size_heatmap = ggplot(beta_df, aes(x = figure_name, y = gene_name, fill = beta_quantile)) + 
    facet_wrap(~phenotype) + 
    geom_tile() + 
    ylab(ylabel) + 
    xlab("Condition") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(space = "Lab", low = "#4575B4", 
                         mid = "#FFFFBF", high = "#E24C36", name = "Normalised effect", midpoint = 0) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank())
  return(effect_size_heatmap)
}


#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#### Link eQTLs to caQTLs based on R2 overlap ####
#Import lead eQTL SNPs
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_qtl_df = purrr::map_df(rasqual_min_pvalues, 
                               ~dplyr::filter(., p_eigen < fdr_thresh), .id = "condition_name") %>% 
  dplyr::arrange(p_nominal)
joint_pairs = dplyr::select(rasqual_qtl_df, gene_id, snp_id) %>% unique() 

#Import ATAC QTL variants
atac_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
atac_qtl_df = purrr::map_df(atac_min_pvalues, 
                            ~dplyr::filter(., p_eigen < fdr_thresh), .id = "condition_name") %>% 
  dplyr::arrange(p_nominal)
atac_joint_pairs = dplyr::select(atac_qtl_df, gene_id, snp_id) %>% unique() 

#Find overlaps using the GWAS overlap code
rna_trait_pairs = addVariantCoords(joint_pairs, vcf_file$snpspos)
atac_trait_pairs = addVariantCoords(atac_joint_pairs, vcf_file$snpspos) %>%
  dplyr::rename(peak_id = gene_id)

#Find overlaps chr by chr
chr_list = idVectorToList(unique(rna_trait_pairs$chr))
overlap_list = purrr::map(chr_list, ~findGWASOverlaps(dplyr::filter(rna_trait_pairs, chr == .) %>% 
                                                        dplyr::select(gene_id, snp_id), 
                                                     dplyr::filter(atac_trait_pairs, chr == .), 
                                                     vcf_file, max_distance = 5e5, min_r2 = 0.8))
rna_atac_overlaps = purrr::map_df(overlap_list, identity)
saveRDS(rna_atac_overlaps, "results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")





#Identify shared QTLs
rna_atac_overlaps = readRDS("results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")

#Filter results with coloc
coloc_overlaps = readRDS("results/ATAC_RNA_overlaps/QTL_overlap_list_coloc.rds")
rna_atac_overlaps = dplyr::semi_join(rna_atac_overlaps, coloc_overlaps, by = c("gene_id", "peak_id"))

#Find minimal p-values for each peaks across conditions
atac_unique_pvalues = purrr::map_df(atac_min_pvalues, identity, .id = "condition_name") %>%
  dplyr::filter(gene_id %in% unique(rna_atac_overlaps$peak_id)) %>%
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(peak_id = gene_id, p_nominal)

#Find unique pairs between genes and peaks
unique_pairs_r2 = dplyr::left_join(rna_atac_overlaps, atac_unique_pvalues, by = "peak_id") %>% 
  dplyr::group_by(gene_id, snp_id) %>% 
  dplyr::arrange(gene_id, snp_id, -R2, p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(chr != "X") %>%
  dplyr::ungroup()

#Fetch summary stats for all QTL pairs
#Import selected p-values from disk
rna_selected_pvalues = readRDS("results/SL1344/eQTLs/rasqual_selected_pvalues.rds")
atac_selected_pvalues = fetchRasqualSNPs(unique_pairs_r2$snp_id, vcf_file$snpspos, qtlResults()$atac_rasqual)

#Import eQTL clusters
variable_qtls = readRDS("results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")
gene_clusters = dplyr::select(variable_qtls$appear, gene_id, snp_id, new_cluster_id) %>% ungroup() %>% unique()

#Extract individual gene clusters
gene_cluster_list = list(IFNg = dplyr::filter(gene_clusters, new_cluster_id %in% c(5,6)),
                         SL1344 = dplyr::filter(gene_clusters, new_cluster_id %in% c(3,4)),
                         IFNg_SL1344 = dplyr::filter(gene_clusters, new_cluster_id %in% c(1)))
gene_cluster_conditions = list(IFNg = c("naive","IFNg"), SL1344 = c("naive","SL1344"),
                               IFNg_SL1344 = c("naive","IFNg", "SL1344", "IFNg_SL1344"))

#Extract cluster pairs
pairs_list = purrr::map(gene_cluster_list, ~extractGenePeakPairs(unique_pairs_r2, .))

#Extract betas for all pairs
betas_list = purrr::map(pairs_list, ~extractBetasForQTLPairs(., rna_selected_pvalues, atac_selected_pvalues) %>%
                          betaCorrectSignPairs())

#Filter and process the betas for plotting
beta_processed = purrr::map2(betas_list, gene_cluster_conditions, ~dplyr::filter(.x, condition_name %in% .y) %>%
                               dplyr::left_join(figureNames(), by = "condition_name") %>%
                               dplyr::mutate(beta_quantile = quantileNormaliseBeta(beta)) %>%
                               dplyr::left_join(gene_name_map, by = "gene_id") %>%
                               sortByBeta("ATAC"))

#Make a heatmaps
plotQTLBetas(beta_processed$IFNg)
plotQTLBetas(beta_processed$SL1344)
plotQTLBetas(beta_processed$IFNg_SL1344)

#Make a line plot
ggplot(beta_processed$IFNg, aes(x = figure_name, y = beta_quantile, group = gene_id)) + 
  geom_point() + geom_line() + facet_wrap(~phenotype)
ggplot(beta_processed$SL1344, aes(x = figure_name, y = beta_quantile, group = gene_id)) + 
  geom_point() + geom_line() + facet_wrap(~phenotype)
ggplot(beta_processed$IFNg_SL1344, aes(x = figure_name, y = beta_quantile, group = gene_id)) + 
  geom_point() + geom_line() + facet_wrap(~phenotype)





###### Reverse analysis #####
#Find minimal p-values for each gene across conditions
rna_unique_pvalues = purrr::map_df(rasqual_min_pvalues, identity, .id = "condition_name") %>%
  dplyr::filter(gene_id %in% unique(rna_atac_overlaps$gene_id)) %>%
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(gene_id, p_nominal)

#Find unique pairs between genes and peaks (focussing on peaks)
atac_unique_pairs_r2 = dplyr::left_join(rna_atac_overlaps, rna_unique_pvalues, by = "gene_id") %>% 
  dplyr::group_by(peak_id, snp_id) %>% 
  dplyr::arrange(peak_id, snp_id, -R2, p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(chr != "X") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(snp_id = gwas_snp_id)

#Import condition-specific QTLs
atac_variable_qtls = readRDS("results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")
peak_clusters = dplyr::select(atac_variable_qtls$appear, gene_id, snp_id, new_cluster_id) %>% ungroup() %>% unique()

#Extract individual gene clusters
peak_cluster_list = list(IFNg = dplyr::filter(peak_clusters, new_cluster_id %in% c(5,6)),
                         SL1344 = dplyr::filter(peak_clusters, new_cluster_id %in% c(3,4)),
                         IFNg_SL1344 = dplyr::filter(peak_clusters, new_cluster_id %in% c(1)))
peak_cluster_conditions = list(IFNg = c("naive","IFNg"), SL1344 = c("naive","SL1344"),
                               IFNg_SL1344 = c("naive","IFNg", "SL1344", "IFNg_SL1344"))


#Import RNA selected p-values for ATAC QTLs
atac_atac_selected_pvalues = readRDS("results/ATAC/QTLs/rasqual_selected_pvalues.rds")
atac_rna_selected_pvalues = fetchRasqualSNPs(atac_unique_pairs_r2$snp_id, vcf_file$snpspos, qtlResults()$rna_rasqual)

#Extract cluster pairs
peak_pairs_list = purrr::map(peak_cluster_list, ~extractPeakGenePairs(atac_unique_pairs_r2, .))

#Extract betas for all pairs
peak_betas_list = purrr::map(peak_pairs_list, ~extractBetasForQTLPairs(., atac_rna_selected_pvalues, atac_atac_selected_pvalues) %>%
                          betaCorrectSignPairs())

#Filter and process the betas for plotting
peak_beta_processed = purrr::map2(peak_betas_list, peak_cluster_conditions, ~dplyr::filter(.x, condition_name %in% .y) %>%
                               dplyr::left_join(figureNames(), by = "condition_name") %>%
                               dplyr::mutate(beta_quantile = quantileNormaliseBeta(beta)) %>%
                               dplyr::mutate(gene_name = peak_id) %>%
                               sortByBeta("RNA"))

#Make a heatmaps
plotQTLBetas(peak_beta_processed$IFNg)
plotQTLBetas(peak_beta_processed$SL1344)
plotQTLBetas(peak_beta_processed$IFNg_SL1344)


ggplot(peak_beta_processed$IFNg, aes(x = figure_name, y = beta_quantile, group = peak_id)) + 
  geom_point() + geom_line() + facet_wrap(~phenotype)
ggplot(peak_beta_processed$SL1344, aes(x = figure_name, y = beta_quantile, group = peak_id)) + 
  geom_point() + geom_line() + facet_wrap(~phenotype)
ggplot(peak_beta_processed$IFNg_SL1344, aes(x = figure_name, y = beta_quantile, group = peak_id)) + 
  geom_point() + geom_line() + facet_wrap(~phenotype)
