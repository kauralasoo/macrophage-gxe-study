library("devtools")
library("dplyr")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("../seqUtils/")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import eQTL minimal p-values
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
fastqtl_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")
fastqtl_100kb = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues_100kb.rds")

#Import cQTL minimal p-values
atac_rasqual_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
atac_fastqtl_min_pvalues = readRDS("results/ATAC/QTLs/fastqtl_min_pvalues.rds")

#Import tQTL results
salmon_min_pvalues = readRDS("results/SL1344/salmon/salmon_qtl_hits.rds")
#Aggregare reviseAnnotations at gene level
salmon_gene_pvalues = purrr::map(salmon_min_pvalues, ~dplyr::arrange(., ensembl_gene_id, p_nominal) %>% 
                                   dplyr::group_by(ensembl_gene_id) %>% 
                                   dplyr::filter(row_number() == 1) %>%
                                   dplyr::ungroup() %>%
                                   dplyr::mutate(gene_id = ensembl_gene_id))
salmon_ensembl_pvalues = readRDS("results/SL1344/salmon/salmon_ensembl_qtl_hits.rds")
leafcutter_min_pvalues = readRDS("results/SL1344/leafcutter/leafcutter_cluster_min_pvalues.rds") %>%
  purrr::map(~dplyr::rename(.,gene_id = cluster_id))
leafcutter_gene_pvalues = purrr::map(leafcutter_min_pvalues, ~dplyr::arrange(., ensembl_gene_id, p_nominal) %>%
                                   dplyr::filter(!is.na(ensembl_gene_id)) %>%
                                   dplyr::group_by(ensembl_gene_id) %>% 
                                   dplyr::filter(row_number() == 1) %>%
                                   dplyr::ungroup() %>%
                                   dplyr::mutate(gene_id = ensembl_gene_id))

#Calculate concordance between conditions
salmon_event_concordance = calculatePairwiseConcordance(salmon_min_pvalues, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "reviseAnnotations", class = "Transcript ratio")
salmon_gene_concordance = calculatePairwiseConcordance(salmon_gene_pvalues, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second)
salmon_ensembl_concordance = calculatePairwiseConcordance(salmon_ensembl_pvalues, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "Ensembl_85", class = "Transcript ratio")
leafcutter_concordance = calculatePairwiseConcordance(leafcutter_min_pvalues, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "LeafCutter", class = "Transcript ratio")
leafcutter_gene_concordance = calculatePairwiseConcordance(leafcutter_gene_pvalues, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second)

#Calculate concordance for eQTLs between conditions
fastqtl_100kb_concordance_matrix = calculatePairwiseConcordance(fastqtl_100kb, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (100kb)", class = "Gene expression")
fastqtl_500kb_concordance_matrix = calculatePairwiseConcordance(fastqtl_min_pvalues, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second)  %>% dplyr::mutate(method = "FastQTL (500kb)", class = "Gene expression")
rasqual_500kb_concordance_matrix = calculatePairwiseConcordance(rasqual_min_pvalues, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second)  %>% dplyr::mutate(method = "RASQUAL (500kb)", class = "Gene expression")

#Calculate concordance for caQTLs between conditions
atac_rasqual_concordance = calculatePairwiseConcordance(atac_rasqual_min_pvalues, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "RASQUAL (50kb)", class = "Chromatin")
atac_fastqtl_concordance = calculatePairwiseConcordance(atac_fastqtl_min_pvalues, ~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE) %>%
  dplyr::filter(first != second) %>% dplyr::mutate(method = "FastQTL (50kb)", class = "Chromatin")


#Make a concordance boxplot for Chapter 3
condition_specificity_df = dplyr::bind_rows(fastqtl_100kb_concordance_matrix, fastqtl_500kb_concordance_matrix, rasqual_500kb_concordance_matrix,
                                            salmon_event_concordance, salmon_ensembl_concordance, leafcutter_concordance)
specificity_plot = ggplot(condition_specificity_df, aes(x = method, y = concordance)) + 
  facet_wrap(~class, scales = "free_x") + 
  geom_boxplot() + geom_point() + theme_light() +
  ylab("Lead variant concordance") +
  scale_y_continuous(limits = c(0,0.8)) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
ggsave("figures/supplementary/concordance_between_conditions.pdf", specificity_plot, width = 6, height = 4)

#Make a concordance boxplot for Chapter 4
ch4_condition_specificity_df = dplyr::bind_rows(fastqtl_500kb_concordance_matrix, rasqual_500kb_concordance_matrix,
                                            atac_rasqual_concordance, atac_fastqtl_concordance) %>%
  dplyr::mutate(class = factor(class, levels = c("Gene expression","Chromatin")))
ch4_specificity_plot = ggplot(ch4_condition_specificity_df, aes(x = method, y = concordance)) + 
  facet_wrap(~class, scales = "free_x") + 
  geom_boxplot() + geom_point() + theme_light() +
  ylab("Lead variant concordance") +
  scale_y_continuous(limits = c(0,0.8)) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
ggsave("figures/supplementary/concordance_between_conditions_atac.pdf", ch4_specificity_plot, width = 6, height = 4)


#Concordance between RASQUAL and FastQTL results
paired_list = purrr::map2(rasqual_min_pvalues, fastqtl_min_pvalues, function(x,y){return(list(RASQUAL = x, FastQTL = y))})
rasqual_fastqtl_concordance = purrr::map(paired_list, 
      ~calculatePairwiseConcordance(.,~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE, R2_thresh = 0.8)) %>%
  purrr::map_df(identity, .id = "condition_name") %>% 
  dplyr::filter(first != second) %>% 
  dplyr::mutate(comparison = paste(first, second, sep = " - "))

#Concordance between tQTLs and eQTLs
paired_list = purrr::map2(fastqtl_100kb, salmon_ensembl_pvalues, function(x,y){return(list(FastQTL = x, Ensembl_85 = y))})
fastqtl_ensembl85_pairwise_concordance = purrr::map(paired_list, 
                                            ~calculatePairwiseConcordance(.,~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE, R2_thresh = 0.8)) %>%
  purrr::map_df(identity, .id = "condition_name") %>% 
  dplyr::filter(first != second) %>% 
  dplyr::mutate(comparison = paste(first, second, sep = " - "))

#Concordance between tQTLs and eQTLs
paired_list = purrr::map2(fastqtl_100kb, salmon_gene_pvalues, function(x,y){return(list(FastQTL = x, reviseAnnotations = y))})
eQTL_tQTL_pairwise_concordance = purrr::map(paired_list, 
      ~calculatePairwiseConcordance(.,~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE, R2_thresh = 0.8)) %>%
  purrr::map_df(identity, .id = "condition_name") %>% 
  dplyr::filter(first != second) %>% 
  dplyr::mutate(comparison = paste(first, second, sep = " - "))

#Concordance between fastqtl and leafCutter estimates
paired_list = purrr::map2(fastqtl_100kb, leafcutter_gene_pvalues, function(x,y){return(list(FastQTL = x, LeafCutter = y))})
FastQTL_leafcutter_concordance = purrr::map(paired_list, 
                                            ~calculatePairwiseConcordance(.,~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE, R2_thresh = 0.8)) %>%
  purrr::map_df(identity, .id = "condition_name") %>% 
  dplyr::filter(first != second) %>% 
  dplyr::mutate(comparison = paste(first, second, sep = " - "))


#Concordance between revised and raw transcript annotations
paired_list = purrr::map2(salmon_gene_pvalues, salmon_ensembl_pvalues, function(x,y){return(list(reviseAnnotations = x, Ensembl_85 = y))})
revised_ensembl_concordance = purrr::map(paired_list, 
      ~calculatePairwiseConcordance(.,~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE, R2_thresh = 0.8)) %>%
  purrr::map_df(identity, .id = "condition_name") %>% 
  dplyr::filter(first != second) %>% 
  dplyr::mutate(comparison = paste(first, second, sep = " - "))

#Concordance between Salmon and leafCutter estimates
paired_list = purrr::map2(salmon_gene_pvalues, leafcutter_gene_pvalues, function(x,y){return(list(reviseAnnotations = x, LeafCutter = y))})
revised_leafcutter_concordance = purrr::map(paired_list, 
                                         ~calculatePairwiseConcordance(.,~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE, R2_thresh = 0.8)) %>%
  purrr::map_df(identity, .id = "condition_name") %>% 
  dplyr::filter(first != second) %>% 
  dplyr::mutate(comparison = paste(first, second, sep = " - "))

#Concordance between RASQUAL and FastQTL results for ATAC-seq
paired_list = purrr::map2(atac_rasqual_min_pvalues, atac_fastqtl_min_pvalues, function(x,y){return(list(RASQUAL = x, FastQTL = y))})
atac_rasqual_fastqtl_concordance = purrr::map(paired_list, 
                                         ~calculatePairwiseConcordance(.,~p_fdr < 0.01, vcf_file$genotypes, tidy = TRUE, R2_thresh = 0.8)) %>%
  purrr::map_df(identity, .id = "condition_name") %>% 
  dplyr::filter(first != second) %>% 
  dplyr::mutate(comparison = paste(first, second, sep = " - "))


#Merge all comparisons together
concordance_df = dplyr::bind_rows(rasqual_fastqtl_concordance, revised_ensembl_concordance, revised_leafcutter_concordance,
                                  eQTL_tQTL_pairwise_concordance, FastQTL_leafcutter_concordance) %>%
  dplyr::mutate(comparison = factor(comparison, levels = unique(comparison))) %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344"))) %>%
  dplyr::mutate(type = c(rep("eQTLs",8), rep("trQTLs",16), rep("eQTLs vs trQTLs",16))) %>%
  dplyr::mutate(type = factor(type, levels = unique(type)))

concordance_plot = ggplot(concordance_df, aes(x = comparison, y = concordance, color = condition_name)) + 
  geom_point() + 
  theme_light() +
  facet_grid(.~ type, scales = "free_x", space = "free_x") +
  scale_color_manual(values = conditionPalette()) + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_y_continuous(limit = c(0,0.7))
ggsave("figures/supplementary/concordance_between_methods.pdf", concordance_plot, width = 8, height = 5)


#Make concordance plot for the caQTL data
concordance_df = dplyr::bind_rows(rasqual_fastqtl_concordance, atac_rasqual_fastqtl_concordance) %>%
  dplyr::mutate(comparison = factor(comparison, levels = unique(comparison))) %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344"))) %>%
  dplyr::mutate(type = c(rep("eQTLs",8), rep("caQTLs",8))) %>%
  dplyr::mutate(type = factor(type, levels = unique(type)))

concordance_plot = ggplot(concordance_df, aes(x = comparison, y = concordance, color = condition_name)) + 
  geom_point(position = position_jitter(width = 0.2)) + 
  theme_light() +
  facet_grid(.~ type, scales = "free_x", space = "free_x") +
  scale_color_manual(values = conditionPalette()) + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_y_continuous(limit = c(0,0.9))
ggsave("figures/supplementary/concordance_between_methods_atac.pdf", concordance_plot, width = 8, height = 5)


