library("devtools")
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
    dplyr::arrange(gene_id, p_nominal) %>%
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
    dplyr::arrange(peak_id, p_nominal) %>%
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

sortByBeta <- function(beta_df, phenotype_name, beta_thresh = 0.59){
  sorted_names = dplyr::filter(beta_df, phenotype == phenotype_name) %>% 
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(beta_mean = mean(beta)) %>%
    dplyr::mutate(beta_diff = beta[2] - beta[1]) %>%
    dplyr::ungroup() %>%
    dplyr::filter(condition_name == "naive") %>%
    dplyr::arrange(beta_quantile) %>%
    dplyr::mutate(type = ifelse(beta > beta_thresh, "Indirect", "Direct"))
  effect_type_df = dplyr::select(sorted_names, peak_id, type, beta_diff) %>% unique()
  sorted_df = dplyr::mutate(beta_df, gene_name = factor(as.character(gene_name), 
                                                        levels = as.character(sorted_names$gene_name))) %>%
    dplyr::left_join(effect_type_df, by = "peak_id")
  return(sorted_df)
}

plotQTLBetasAll2 <- function(beta_df){
  n_pairs = nrow(dplyr::select(beta_df, gene_name, snp_id) %>% unique())
  ylabel = paste(n_pairs, "eQTL-caQTL pairs")
  effect_size_heatmap = ggplot(beta_df, aes(x = phenotype, y = qtl_id, fill = beta_quantile)) + 
    facet_wrap(~figure_name) + 
    geom_tile() + 
    ylab(ylabel) + 
    xlab("Condition") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(space = "Lab", low = "#4575B4", 
                         mid = "#FFFFBF", high = "#E24C36", name = "Normalised effect", midpoint = 0) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank()) +
    theme(legend.title = element_text(angle = 90)) + 
    theme(axis.text.x = element_text(angle = 15))
  return(effect_size_heatmap)
}

forwardAnalysis <- function(atac_min_pvalues, rna_atac_overlaps, use_filtering = FALSE, filter_threshold = 1){
  
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
    dplyr::arrange(gene_id, snp_id, p_nominal) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::filter(chr != "X") %>%
    dplyr::ungroup()
  
  #Fetch summary stats for all QTL pairs
  #Import selected p-values from disk
  rna_selected_pvalues = readRDS("results/SL1344/eQTLs/rasqual_selected_pvalues.rds")
  atac_selected_pvalues = fetchRasqualSNPs(unique_pairs_r2$snp_id, vcf_file$snpspos, qtlResults()$atac_rasqual)
  
  #Import eQTL clusters
  variable_qtls = readRDS("results/SL1344/eQTLs/appear_disappear_eQTLs.rds")
  gene_clusters = dplyr::select(variable_qtls$appear, gene_id, snp_id, max_condition) %>% ungroup() %>% unique()
  
  #Perform additional, more stringent filtering on the eQTL effect size
  if(use_filtering == TRUE){
    gene_cluster_effects = dplyr::filter(variable_qtls$appear, condition_name == "naive" | condition_name == max_condition) %>% 
      dplyr::transmute(gene_id, snp_id, max_condition, condition_name, beta) %>%
      dplyr::mutate(condition_name = ifelse(condition_name == "naive", "naive", "max")) %>%
      tidyr::spread(condition_name, beta)
    gene_clusters = dplyr::filter(gene_cluster_effects, abs(max) >= filter_threshold, 
                                  abs(naive) <= filter_threshold, 
                                  abs(max-naive) >= filter_threshold) %>%
      dplyr::select(gene_id, snp_id, max_condition)
  }
  
  #Extract individual gene clusters
  gene_cluster_list = list(IFNg = dplyr::filter(gene_clusters, max_condition == "IFNg"),
                           SL1344 = dplyr::filter(gene_clusters, max_condition == "SL1344"),
                           IFNg_SL1344 = dplyr::filter(gene_clusters, max_condition == "IFNg_SL1344"))
  gene_cluster_conditions = list(IFNg = c("naive","IFNg"), SL1344 = c("naive","SL1344"),
                                 IFNg_SL1344 = c("naive","IFNg_SL1344"))
  
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
                                 sortByBeta("ATAC") %>%
                                 dplyr::mutate(phenotype = ifelse(phenotype == "ATAC", "ATAC-seq", "RNA-seq")))
  
  #Merge all betas together
  all_betas = purrr::map_df(beta_processed, ~dplyr::arrange(., gene_name) %>%
                              dplyr::mutate(.,stimulation_state = ifelse(condition_name == "naive", "Naive","Stimulated")), .id = "max_effect") %>%
    dplyr::mutate(qtl_id = paste(gene_name, snp_id)) %>% 
    dplyr::mutate(qtl_id = factor(qtl_id, levels = unique(qtl_id))) %>%
    dplyr::mutate(max_effect = ifelse(max_effect == "IFNg", "I", ifelse(max_effect == "SL1344", "S", "I+S"))) %>%
    dplyr::mutate(max_effect = factor(max_effect, levels = c("I","S","I+S")))
  
  #Count caQTLs present in the naive condition
  present_fraction = dplyr::filter(all_betas, phenotype == "ATAC-seq", condition_name == "naive") %>% 
    dplyr::mutate(beta_binary = ifelse(abs(beta) > filter_threshold, "present", "absent")) %>% 
    dplyr::group_by(max_effect, beta_binary) %>% 
    dplyr::summarise(count = length(beta_binary)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(beta_binary, count) %>%
    dplyr::mutate(fraction = present/(absent+present)) %>%
    dplyr::mutate(type = "forward")
  
  return(list(betas = all_betas, present_fraction = present_fraction))
}

reverseAnalysis <- function(rna_min_pvalues, rna_atac_overlaps, use_filtering = FALSE, filter_threshold = 1){
  ###### Reverse analysis #####
  #Find minimal p-values for each gene across conditions
  rna_unique_pvalues = purrr::map_df(rna_min_pvalues, identity, .id = "condition_name") %>%
    dplyr::filter(gene_id %in% unique(rna_atac_overlaps$gene_id)) %>%
    dplyr::group_by(gene_id) %>% 
    dplyr::arrange(gene_id, p_nominal) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(gene_id, p_nominal)
  
  #Find unique pairs between genes and peaks (focussing on peaks)
  atac_unique_pairs_r2 = dplyr::left_join(rna_atac_overlaps, rna_unique_pvalues, by = "gene_id") %>% 
    dplyr::group_by(peak_id, snp_id) %>% 
    dplyr::arrange(peak_id, snp_id, p_nominal) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::filter(chr != "X") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(snp_id = gwas_snp_id)
  
  #Import condition-specific QTLs
  atac_variable_qtls = readRDS("results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")
  peak_clusters = dplyr::select(atac_variable_qtls$appear, gene_id, snp_id, max_condition) %>% ungroup() %>% unique()
  
  #Perform additional, more stringent filtering on the eQTL effect size
  if(use_filtering == TRUE){
    peak_cluster_effects = dplyr::filter(atac_variable_qtls$appear, condition_name == "naive" | condition_name == max_condition) %>% 
      dplyr::transmute(gene_id, snp_id, max_condition, condition_name, beta) %>%
      dplyr::mutate(condition_name = ifelse(condition_name == "naive", "naive", "max")) %>%
      tidyr::spread(condition_name, beta)
    peak_clusters = dplyr::filter(peak_cluster_effects, abs(max) >= filter_threshold, 
                                  abs(naive) <= filter_threshold, 
                                  abs(max-naive) >= filter_threshold) %>%
      dplyr::select(gene_id, snp_id, max_condition)
  }
  
  #Extract individual gene clusters
  peak_cluster_list = list(IFNg = dplyr::filter(peak_clusters, max_condition == "IFNg"),
                           SL1344 = dplyr::filter(peak_clusters, max_condition == "SL1344"),
                           IFNg_SL1344 = dplyr::filter(peak_clusters, max_condition == "IFNg_SL1344"))
  peak_cluster_conditions = list(IFNg = c("naive","IFNg"), SL1344 = c("naive","SL1344"),
                                 IFNg_SL1344 = c("naive","IFNg_SL1344"))
  
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
                                      sortByBeta("RNA") %>% 
                                      dplyr::mutate(phenotype = ifelse(phenotype == "ATAC", "ATAC-seq", "RNA-seq")))
  
  #Merge all betas together
  peak_all_betas = purrr::map_df(peak_beta_processed, ~dplyr::arrange(., gene_name) %>%
                                   dplyr::mutate(.,stimulation_state = ifelse(condition_name == "naive", "Naive","Stimulated")), .id = "max_effect") %>%
    dplyr::mutate(qtl_id = paste(gene_name, snp_id)) %>% 
    dplyr::mutate(qtl_id = factor(qtl_id, levels = unique(qtl_id))) %>%
    dplyr::mutate(max_effect = ifelse(max_effect == "IFNg", "I", ifelse(max_effect == "SL1344", "S", "I+S"))) %>%
    dplyr::mutate(max_effect = factor(max_effect, levels = c("I","S","I+S")))
  
  #Count caQTLs present in the naive condition
  peak_present_fraction = dplyr::filter(peak_all_betas, phenotype == "RNA-seq", condition_name == "naive") %>% 
    dplyr::mutate(beta_binary = ifelse(abs(beta) > filter_threshold, "present", "absent")) %>% 
    dplyr::group_by(max_effect, beta_binary) %>% 
    dplyr::summarise(count = length(beta_binary)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(beta_binary, count) %>%
    dplyr::mutate(fraction = present/(absent+present)) %>%
    dplyr::mutate(type = "reverse")
  
  return(list(all_betas = peak_all_betas, present_fraction = peak_present_fraction))
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

#Import shared QTLs
rna_atac_overlaps = readRDS("results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")



#Original analysis (FC threshold 0.59)
#Perform forward and reverse analysis indepedndently
fwd_results = forwardAnalysis(atac_min_pvalues, rna_atac_overlaps, use_filtering = FALSE, filter_threshold = 0.59)
rev_results = reverseAnalysis(rasqual_min_pvalues, rna_atac_overlaps, use_filtering = FALSE, filter_threshold = 0.59)

#Make a plot estimating the proportion of foreshadowing
combined_results = dplyr::bind_rows(fwd_results$present_fraction, rev_results$present_fraction) %>%
  dplyr::mutate(present = ifelse(is.na(present),0,present)) %>%
  dplyr::mutate(fraction = ifelse(is.na(fraction),0,fraction))

plot_data = combined_results %>%
  dplyr::mutate(type = ifelse(type == "forward", "caQTL before eQTL", "eQTL before caQTL"))

foreshadow_plot = ggplot(plot_data, aes(x = max_effect, y = fraction, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  xlab("Condition") +
  ylab("Fraction of caQTL-eQTL pairs") + 
  theme_light() +
  theme(legend.position = "right")

#Save results to disk
ggsave("figures/main_figures/foreshadowing_proportions.pdf", plot = foreshadow_plot, width = 3.7, height = 3)
write.table(combined_results, "results/ATAC_RNA_overlaps/foreshadow_quant.txt", sep = "\t", quote = FALSE)


#Analysis with a more stringent cutoff (FC threshold 1)
#Perform forward and reverse analysis indepedndently
fwd_results = forwardAnalysis(atac_min_pvalues, rna_atac_overlaps, use_filtering = TRUE, filter_threshold = 1)
rev_results = reverseAnalysis(rasqual_min_pvalues, rna_atac_overlaps, use_filtering = TRUE, filter_threshold = 1)

#Make a plot estimating the proportion of foreshadowing
combined_results = dplyr::bind_rows(fwd_results$present_fraction, rev_results$present_fraction) %>%
  dplyr::mutate(present = ifelse(is.na(present),0,present)) %>%
  dplyr::mutate(fraction = ifelse(is.na(fraction),0,fraction))

plot_data = combined_results %>%
  dplyr::mutate(type = ifelse(type == "forward", "caQTL before eQTL", "eQTL before caQTL"))

foreshadow_plot = ggplot(plot_data, aes(x = max_effect, y = fraction, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  xlab("Condition") +
  ylab("Fraction of caQTL-eQTL pairs") + 
  theme_light() +
  theme(legend.position = "right")

ggsave("figures/supplementary/foreshadowing_proportions.FC_2.pdf", plot = foreshadow_plot, width = 3.7, height = 3)
ggsave("figures/supplementary/foreshadowing_proportions.FC_2.png", plot = foreshadow_plot, width = 3.7, height = 3)


#Perform the same analysis using coloc pairs

#Filter results with coloc
coloc_overlaps = readRDS("results/ATAC_RNA_overlaps/QTL_overlap_list_coloc.rds")
rna_atac_overlaps = dplyr::semi_join(rna_atac_overlaps, coloc_overlaps, by = c("gene_id", "peak_id"))

#Perform forward and reverse analysis
fwd_results = forwardAnalysis(atac_min_pvalues, rna_atac_overlaps, use_filtering = FALSE, filter_threshold = 0.59)
rev_results = reverseAnalysis(rasqual_min_pvalues, rna_atac_overlaps, use_filtering = FALSE, filter_threshold = 0.59)

#Count fraction present
combined_results = dplyr::bind_rows(fwd_results$present_fraction, rev_results$present_fraction) %>%
  dplyr::mutate(present = ifelse(is.na(present),0,present)) %>%
  dplyr::mutate(absent = ifelse(is.na(absent),0,absent)) %>%
  dplyr::mutate(fraction = present/(absent+present))

joint_counts = combined_results %>% 
  dplyr::group_by(type) %>% 
  dplyr::summarise(absent_sum = sum(absent), present_sum = sum(present)) %>%
  dplyr::mutate(fraction = present_sum/(present_sum + absent_sum))

plot_data = joint_counts %>%
  dplyr::mutate(type = ifelse(type == "forward", "caQTL before eQTL", "eQTL before caQTL")) %>%
  dplyr::mutate(max_effect = "All conditions")

foreshadow_plot = ggplot(plot_data, aes(x = max_effect, y = fraction, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  xlab("Condition") +
  ylab("Fraction of caQTL-eQTL pairs") + 
  theme_light() +
  theme(legend.position = "right")

#Save files
ggsave("figures/supplementary/foreshadowing_proportions.coloc.pdf", plot = foreshadow_plot, width = 3.7, height = 3)
write.table(joint_counts, "figures/tables/foreshadow_quant.coloc.txt", sep = "\t", quote = FALSE)



#Perform downsampling analyis

#Perform down-sampling anaysis on the QTLs
performSubsamplingAnalysis <- function(rna_atac_overlaps, rasqual_qtl_df, atac_qtl_df, atac_min_pvalues, rasqual_min_pvalues){
  print("item")
  uniq_genes = unique(rasqual_qtl_df$gene_id)
  uniq_peaks = unique(atac_qtl_df$gene_id)
  sample_peaks = sample(uniq_peaks, length(uniq_genes))
  rna_atac_overlaps_sample = dplyr::filter(rna_atac_overlaps, peak_id %in% sample_peaks)
  
  #Perform FWD and REV analysis
  fwd_results = forwardAnalysis(atac_min_pvalues, rna_atac_overlaps_sample, use_filtering = FALSE, filter_threshold = 0.59)
  rev_results = reverseAnalysis(rasqual_min_pvalues, rna_atac_overlaps_sample, use_filtering = FALSE, filter_threshold = 0.59)
  
  #Make a plot estimating the proportion of foreshadowing
  joint_counts = dplyr::bind_rows(fwd_results$present_fraction, rev_results$present_fraction) %>%
    dplyr::mutate(present = ifelse(is.na(present),0,present)) %>%
    dplyr::mutate(absent = ifelse(is.na(absent),0,absent)) %>%
    dplyr::mutate(fraction = present/(absent+present)) %>%
    dplyr::group_by(type) %>% 
    dplyr::summarise(absent_sum = sum(absent), present_sum = sum(present)) %>%
    dplyr::mutate(fraction = present_sum/(present_sum + absent_sum))
  
  return(joint_counts)
}

#Apply to 100 random samples
list = as.list(c(1:100))
quietSampling = purrr::quietly(performSubsamplingAnalysis)
result = purrr::map(list, ~quietSampling(rna_atac_overlaps, rasqual_qtl_df, 
                                                      atac_qtl_df, atac_min_pvalues, rasqual_min_pvalues))
df_list = purrr::map(result, ~.$result)
names(df_list) = as.character(c(1:100))

fraction_df = purrr::map_df(df_list, identity, .id = "sample") %>%
  dplyr::mutate(type = ifelse(type == "forward", "caQTL before eQTL", "eQTL before caQTL"))
plot = ggplot(fraction_df, aes(x = fraction, fill = type)) + 
  geom_histogram(binwidth = .03) +
  theme_light() +
  xlab("Fraction of caQTL-eQTL pairs")
ggsave("figures/supplementary/foreshadowing_proportions.downsample_caQTLs.pdf", plot = plot, width = 4, height = 3)
ggsave("figures/supplementary/foreshadowing_proportions.downsample_caQTLs.png", plot = plot, width = 4, height = 3)




all_betas_plot1 = plotQTLBetasAll2(dplyr::filter(all_betas, max_effect %in% c("I")))
all_betas_plot2 = plotQTLBetasAll2(dplyr::filter(all_betas, max_effect %in% c("S")))
all_betas_plot3 = plotQTLBetasAll2(dplyr::filter(all_betas, max_effect %in% c("I+S")))

dplyr::filter(all_betas, max_effect %in% c("I+S")) %>% 
  dplyr::filter(phenotype == "ATAC-seq", condition_name == "naive") %>% 
  dplyr::arrange(-abs(beta)) %>% 
  dplyr::filter(abs(beta) > 0.59) %>% 
  dplyr::select(gene_name, beta) %>%
  tail()

#Plot with gene names
all_betas_plot1 + theme(axis.text.y = element_text())
all_betas_plot2 + theme(axis.text.y = element_text())
all_betas_plot3 + theme(axis.text.y = element_text())

if(use_coloc == TRUE){
  ggsave("figures/supplementary/eQTLs_vs_caQTL_heatmap_1_coloc.pdf", all_betas_plot1, width = 3, height = 2)
  ggsave("figures/supplementary/eQTLs_vs_caQTL_heatmap_2_coloc.pdf", all_betas_plot2, width = 3, height = 2)
  ggsave("figures/supplementary/eQTLs_vs_caQTL_heatmap_3_coloc.pdf", all_betas_plot3, width = 3, height = 4.5)
} else {
  ggsave("figures/main_figures/eQTLs_vs_caQTL_heatmap_1.pdf", all_betas_plot1, width = 3, height = 2)
  ggsave("figures/main_figures/eQTLs_vs_caQTL_heatmap_2.pdf", all_betas_plot2, width = 3, height = 2)
  ggsave("figures/main_figures/eQTLs_vs_caQTL_heatmap_3.pdf", all_betas_plot3, width = 3, height = 4.5)
}




if(use_coloc == FALSE){
  saveRDS(all_betas, "results/ATAC_RNA_overlaps/caQTL_eQTL_pairs_betas.rds")
}





peak_all_betas_plot = plotQTLBetasAll2(peak_all_betas)

if(use_coloc == TRUE){
  ggsave("figures/supplementary/eQTLs_vs_caQTL_heatmap_reverse.coloc.pdf", peak_all_betas_plot, width = 3, height = 7)
} else{
  ggsave("figures/supplementary/eQTLs_vs_caQTL_heatmap_reverse.pdf", peak_all_betas_plot, width = 3, height = 7)
}


if(use_filtering == TRUE){
  plot_data = combined_results %>%
    dplyr::mutate(type = ifelse(type == "forward", "caQTL before eQTL", "eQTL before caQTL"))
  
  foreshadow_plot = ggplot(plot_data, aes(x = max_effect, y = fraction, fill = type)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    xlab("Condition") +
    ylab("Fraction of caQTL-eQTL pairs") + 
    theme_light() +
    theme(legend.position = "right")
  ggsave("figures/supplementary/foreshadowing_proportions.FC_2.pdf", plot = foreshadow_plot, width = 3.7, height = 3)
  ggsave("figures/supplementary/foreshadowing_proportions.FC_2.png", plot = foreshadow_plot, width = 3.7, height = 3)
}


#Count the number of caQTLs linked to eQTLs and vice versa
rna_atac_overlaps = readRDS("results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")
peaks_per_gene = dplyr::select(rna_atac_overlaps, gene_id, peak_id) %>% 
  unique() %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::summarise(peak_count = length(peak_id)) %>%
  dplyr::rename(feature_id = gene_id) %>%
  dplyr::mutate(feature_type = "Number of caQTLs per gene")

genes_per_peak = dplyr::select(rna_atac_overlaps, gene_id, peak_id) %>% 
  unique() %>% 
  dplyr::group_by(peak_id) %>% 
  dplyr::summarise(peak_count = length(gene_id)) %>%
  dplyr::rename(feature_id = peak_id) %>%
  dplyr::mutate(feature_type = "Number of eQTLs per region")

pair_dist = dplyr::bind_rows(peaks_per_gene, genes_per_peak)
pair_plot = ggplot(pair_dist, aes(x = peak_count)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~feature_type, scales = "free_y") +
  scale_x_continuous(limits = c(-0,20)) +
  theme_light() +
  xlab("Number of caQTLs/eQTLs")

ggsave("figures/supplementary/caQTL_eQTL_pairs_histogram.pdf", plot = pair_plot, width = 5, height = 3)
ggsave("figures/supplementary/caQTL_eQTL_pairs_histogram.png", plot = pair_plot, width = 5, height = 3)

#Are master caQTLs more likely to regulate multiple genes?

#Import multi-peak QTLs
result_list = readRDS("results/ATAC/QTLs/qtl_peak_type_assignment.rds")
master_peaks = unique(result_list$dependents$unique_masters$master_id)
master_eQTL_count = dplyr::filter(genes_per_peak, feature_id %in% master_peaks)

fisher.test(matrix(c(4194-3039,3039, 18, 176),ncol = 2 ))



