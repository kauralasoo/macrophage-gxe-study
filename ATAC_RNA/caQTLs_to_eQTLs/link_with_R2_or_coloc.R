library("devtools")
library("plyr")
library("dplyr")
library("purrr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

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

#Construct GRanges object of SNP positions
selected_snps = dplyr::filter(vcf_file$snpspos, snpid %in% unique_pairs_r2$snp_id) %>%
  dplyr::transmute(snp_id = snpid, seqnames = chr, start = pos, end = pos, strand = "*") %>%
  dataFrameToGRanges()

#Fetch corresponding SNPs from ATAC data
atac_tabix_list = qtlResults()$atac_rasqual
atac_selected_pvalues = lapply(atac_tabix_list, function(tabix, snps) rasqualTools::tabixFetchSNPs(snps, tabix), selected_snps)


#Import eQTL clusters
variable_qtls = readRDS("results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")
gene_clusters = dplyr::select(variable_qtls$appear, gene_id, snp_id, new_cluster_id) %>% ungroup() %>% unique()

#Extract IFNg pairs
ifng_genes = dplyr::filter(gene_clusters, new_cluster_id %in% c(5,6))
ifng_pairs = dplyr::semi_join(unique_pairs_r2, ifng_genes, by = c("gene_id", "snp_id")) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(gene_id, -R2, p_nominal) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
ifng_pairs_selected = dplyr::select(ifng_pairs, gene_id, snp_id, peak_id) %>% unique()

#Extract joint betas
joint_betas = extractBetasForQTLPairs(ifng_pairs_selected, rna_selected_pvalues, atac_selected_pvalues) %>%
  betaCorrectSignPairs()

#Make a plot for IFNg pairs
ifng_effects = joint_betas %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  dplyr::group_by(snp_id, peak_id, gene_id) %>% 
  dplyr::mutate(beta_std = (beta - mean(beta))/sd(beta), beta_scaled = beta/max(beta)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_binary = ifelse(beta >= 0.59, 1, 0)) %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg"))) %>%
  dplyr::left_join(gene_name_map, by = "gene_id")

#Caclulate scaled diff
scaled_diff = dplyr::filter(ifng_effects, phenotype == "ATAC") %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(condition_name) %>% 
  dplyr::mutate(scaled_diff = beta_scaled[2] - beta_scaled[1]) %>% 
  dplyr::mutate(mean_beta = mean(beta)) %>% 
  dplyr::filter(condition_name == "naive") %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(beta_scaled)
ifng_effects_sorted = dplyr::mutate(ifng_effects, gene_name = factor(as.character(gene_name), 
                                                                     levels = as.character(scaled_diff$gene_name))) %>%
  dplyr::mutate(phenotype = ifelse(phenotype == "ATAC", "ATAC-seq", "RNA-seq")) %>%
  dplyr::left_join(figureNames(), by = "condition_name")

#Make a heatmap
n_pairs = nrow(dplyr::select(ifng_effects_sorted, gene_name, snp_id) %>% unique())
ylabel = paste(n_pairs, "eQTL-caQTL pairs")
ifng_effect_size_heatmap = ggplot(ifng_effects_sorted, aes(x = figure_name, y = gene_name, fill = beta)) + 
  facet_wrap(~phenotype) + 
  geom_tile() + 
  ylab(ylabel) + 
  xlab("Condition") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", 
                       mid = "#FFFFBF", high = "#E24C36", name = "Relative effect", midpoint = 0) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank())

ggplot(ifng_effects_sorted, aes(x = figure_name, y = beta, group = gene_id)) + 
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
variable_qtls = readRDS("results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")
peak_clusters = dplyr::select(variable_qtls$appear, gene_id, snp_id, new_cluster_id) %>% ungroup() %>% unique()

#Import ATAC selected p-values
atac_selected_pvalues = readRDS("results/ATAC/QTLs/rasqual_selected_pvalues.rds")

#Construct GRanges object of SNP positions
selected_snps = dplyr::filter(vcf_file$snpspos, snpid %in% atac_unique_pairs_r2$snp_id) %>%
                                dplyr::transmute(snp_id = snpid, seqnames = chr, start = pos, end = pos, strand = "*") %>%
                                dataFrameToGRanges()
                                
#Fetch corresponding SNPs from ATAC data
rna_tabix_list = qtlResults()$rna_rasqual
rna_selected_pvalues = lapply(rna_tabix_list, function(tabix, snps) rasqualTools::tabixFetchSNPs(snps, tabix), selected_snps)

#Extract IFNg pairs
ifng_peaks = dplyr::filter(peak_clusters, new_cluster_id %in% c(5,6)) %>%
  dplyr::rename(peak_id = gene_id)
ifng_pairs = dplyr::semi_join(atac_unique_pairs_r2, ifng_peaks, by = c("peak_id", "snp_id")) %>%
  dplyr::group_by(peak_id) %>%
  dplyr::arrange(peak_id, -R2, p_nominal) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
ifng_pairs_selected = dplyr::select(ifng_pairs, gene_id, snp_id, peak_id) %>% unique()

#Extract joint betas
joint_betas = extractBetasForQTLPairs(ifng_pairs_selected, rna_selected_pvalues, atac_selected_pvalues) %>%
  betaCorrectSignPairs()

#Make a plot for IFNg pairs
ifng_effects = joint_betas %>%
  dplyr::filter(condition_name %in% c("naive","IFNg")) %>%
  dplyr::group_by(snp_id, peak_id, gene_id, phenotype) %>% 
  dplyr::mutate(beta_std = (beta - mean(beta))/sd(beta), beta_scaled = beta/max(beta)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_binary = ifelse(beta >= 0.59, 1, 0)) %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg"))) %>%
  dplyr::left_join(gene_name_map, by = "gene_id") 

#Caclulate scaled diff
scaled_diff = dplyr::filter(ifng_effects, phenotype == "RNA") %>% 
  dplyr::group_by(peak_id, gene_id) %>% 
  dplyr::arrange(peak_id,gene_id, condition_name) %>% 
  dplyr::mutate(scaled_diff = beta_scaled[2] - beta_scaled[1]) %>% 
  dplyr::filter(condition_name == "naive") %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(beta)
ifng_effects_sorted = dplyr::mutate(ifng_effects, gene_name = factor(as.character(peak_id), 
                                                                     levels = as.character(scaled_diff$peak_id))) %>%
  dplyr::mutate(phenotype = ifelse(phenotype == "ATAC", "ATAC-seq", "RNA-seq")) %>%
  dplyr::left_join(figureNames(), by = "condition_name")


n_pairs = nrow(dplyr::select(ifng_effects_sorted, gene_name, snp_id) %>% unique())
ylabel = paste(n_pairs, "eQTL-caQTL pairs")
ifng_effect_size_heatmap = ggplot(ifng_effects_sorted, aes(x = figure_name, y = gene_name, fill = beta)) + 
  facet_wrap(~phenotype) + 
  geom_tile() + 
  ylab(ylabel) + 
  xlab("Condition") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", 
                       mid = "#FFFFBF", high = "#D73027", name = "Relative effect", midpoint = 0) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank())

ggplot(ifng_effects_sorted, aes(x = figure_name, y = beta, group = peak_id)) + geom_point() + geom_line() + facet_wrap(~phenotype)

                                
