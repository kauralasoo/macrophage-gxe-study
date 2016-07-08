library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
library("purrr")
library("coloc")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import variant information
snp_info = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

#Load caQTL p-values from disk
atac_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
atac_qtl_df = extractQTLsFromList(atac_min_pvalues, fdr_cutoff = 0.1)

#Load eQTL p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_qtl_df = extractQTLsFromList(rasqual_min_pvalues, fdr_cutoff = 0.1)
joint_pairs = dplyr::select(rasqual_qtl_df, gene_id, snp_id) %>% unique() 
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

#Find unique eQTLs for this gene:
independent_qtls = dplyr::filter(filtered_pairs, gene_id == "ENSG00000120899")
independent_qtl_res = dplyr::semi_join(rasqual_qtl_df, independent_qtls)

#Mark the likely causal variants
causal_variants = c("rs28834970","rs10086852", "rs34181358", "rs56090771")
names(causal_variants) = c("naive", "IFNg_SL1344", "IFNg", "naive2")

#Extract relevant genotypes from vcf
extracted_genotypes = t(vcf_file$genotypes[causal_variants,]) %>% data.frame()
extracted_genotypes = dplyr::mutate(extracted_genotypes, genotype_id = rownames(extracted_genotypes)) %>%
  dplyr::select(genotype_id, everything())

covariate_names = c("PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6", "sex_binary")

#Extract expression data for the gene
qtl_data = constructQtlPlotDataFrame("ENSG00000120899", "rs28834970", combined_expression_data$cqn, vcf_file$genotypes, 
                                       combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::left_join(extracted_genotypes, by = "genotype_id")
naive_plot = plotQtlRow(qtl_data)
ggsave("figures/main_figures/PTK2B_naive_eQTL.pdf", naive_plot, width = 8, height = 4)

ifng_SL1344_plot = constructQtlPlotDataFrame("ENSG00000120899", "rs10086852", combined_expression_data$cqn, vcf_file$genotypes, 
                                     combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::left_join(extracted_genotypes, by = "genotype_id") %>%
  plotQtlRow()
ggsave("figures/main_figures/PTK2B_IFNg_SL1344_eQTL.pdf", ifng_SL1344_plot, width = 8, height = 4)


#Calculate residuals
#Without regerssing out the second QTL
formula = as.formula(paste("norm_exp ~ ",
                           paste(covariate_names, collapse = " + "), sep = "+ "))

data = dplyr::filter(qtl_data, condition_name == "IFNg")
lmfit = lm(formula, data)
data = dplyr::mutate(data, resid = as.numeric(lmfit$residuals))

orig_plot = ggplot2::ggplot(data, ggplot2::aes(x = as.factor(rs28834970) , y = resid , color = condition_name)) + 
  #ggplot2::facet_grid(~condition_name) + 
  ggplot2::geom_boxplot(outlier.shape = NA) + 
  ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
  ggplot2::ylab("Expression residual") +
  ggplot2::theme_light() + 
  ggplot2::scale_color_manual(values = conditionPalette(), guide=FALSE)

summary(lm(resid ~ rs28834970, data))

#With regressing out the second QTL
formula = as.formula(paste("norm_exp ~ rs10086852 ",
                                         paste(covariate_names, collapse = " + "), sep = "+ "))

data = dplyr::filter(qtl_data, condition_name == "IFNg")
lmfit = lm(formula, data)
data = dplyr::mutate(data, resid = as.numeric(lmfit$residuals))

regressed_plot = ggplot2::ggplot(data, ggplot2::aes(x = as.factor(rs28834970) , y = resid , color = condition_name)) + 
  #ggplot2::facet_grid(~condition_name) + 
  ggplot2::geom_boxplot(outlier.shape = NA) + 
  ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
  ggplot2::ylab("Expression residual") +
  ggplot2::theme_light() + 
  ggplot2::scale_color_manual(values = conditionPalette(), guide=FALSE)

summary(lm(resid ~ rs28834970, data))

#Joint the two plots together
joint = cowplot::plot_grid(orig_plot, regressed_plot, labels = c("A","B"))
ggsave("figures/main_figures/PTK2B_sign_flip_residuals.pdf", plot = joint, width = 5, height = 4)



#Fetch all kinds of p-values from the PTK2B region for plotting and coloc
#Define region
ptk2b_region = constructGeneRanges(data_frame(gene_id = "ENSG00000120899"), combined_expression_data$gene_metadata, 5e4)

#Import p-values from the the original GWAS
alzheimer_pvals = importGWASSummaryStats(ptk2b_region, "databases/GWAS/IGAP_summary_statistics/IGAP_stage_1.GRCh38.sorted.bed.gz")[[1]] %>%
  dplyr::select(-snp_id)
alzheimer_pvalues = dplyr::select(snp_info, chr, pos, snp_id, MAF) %>% 
  dplyr::left_join(alzheimer_pvals, ., by = c("chr","pos")) %>% 
  dplyr::filter(!is.na(snp_id)) %>% 
  dplyr::arrange(p_nominal) %>%
  dplyr::mutate(R2 = calculateR2FromLead(snp_id, vcf_file$genotypes)) %>%
  dplyr::mutate(condition_name = "AZ GWAS")

#Import RASQUAL and fastqtl pvalues for PTK2B gene
rna_rasqual_pvalues = purrr::map_df(qtlResults()$rna_rasqual, ~tabixFetchGenes(ptk2b_region, .)[[1]], .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
rna_fastqtl_pvalues = purrr::map_df(qtlResults()$rna_fastqtl, ~fastqtlTabixFetchGenes(ptk2b_region, .)[[1]], .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))

#Add R2 values per conditon
rasqual_r2 = dplyr::group_by(rna_rasqual_pvalues, condition_name) %>% dplyr::arrange(p_nominal) %>%
  dplyr::mutate(R2 = calculateR2FromLead(snp_id, vcf_file$genotypes)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(chr = as.character(chr))

#Merge GWAS and eQTL data together
joint_pvalues = dplyr::bind_rows(dplyr::select(alzheimer_pvalues, chr, pos, snp_id, p_nominal, R2, condition_name),
                 dplyr::select(rasqual_r2, chr, pos, snp_id, p_nominal, R2, condition_name)) %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("AZ GWAS","naive","IFNg","SL1344","IFNg_SL1344")))

#Make a manhattan plot
rasqual_manhattan = ggplot(joint_pvalues, aes(x = pos, y = -log(p_nominal, 10), colour = R2)) + 
  geom_point() + 
  facet_grid(condition_name~., scales = "free_y") +
  theme_light()
ggsave("figures/main_figures/PTK2B_RASQUAL_manhattan.pdf", plot = rasqual_manhattan, width = 6, height = 7)

#Add R2 values per conditon
fastqtl_r2 = dplyr::group_by(rna_fastqtl_pvalues, condition_name) %>% dplyr::arrange(p_nominal) %>%
  dplyr::mutate(R2 = calculateR2FromLead(snp_id, vcf_file$genotypes)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(chr = as.character(chr))

#Merge GWAS and eQTL data together
fastqtl_joint_pvalues = dplyr::bind_rows(dplyr::select(alzheimer_pvalues, chr, pos, snp_id, p_nominal, R2, condition_name),
                                 dplyr::select(fastqtl_r2, chr, pos, snp_id, p_nominal, R2, condition_name)) %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("AZ GWAS","naive","IFNg","SL1344","IFNg_SL1344")))

#Make a manhattan plot
fastqtl_manhattan = ggplot(fastqtl_joint_pvalues, aes(x = pos, y = -log(p_nominal, 10), colour = R2)) + 
  geom_point() + 
  facet_grid(condition_name~., scales = "free_y") +
  theme_light()
ggsave("figures/main_figures/PTK2B_FASTQTL_manhattan.pdf", plot = fastqtl_manhattan, width = 6, height = 7)


#Fetch associatated ATAC peaks for each lead eQTL SNP
peak_midpoints = dplyr::transmute(atac_list$gene_metadata, gene_id, midpoint = end -((end-start)/2)) %>% tbl_df()
atac_snp_results = purrr::map(qtlResults()$atac_rasqual, ~tabixFetchSNPsQuick(causal_variants,.,vcf_file$snpspos)) %>%
  ldply(.id = "condition_name") %>%
  dplyr::left_join(peak_midpoints) %>%
  dplyr::tbl_df()

dplyr::filter(atac_snp_results, p_nominal < 1e-3) %>% dplyr::select(gene_id, snp_id, p_nominal, condition_name)

#Fetch RASQUAL results for the ATAC peak
peak1_region = constructGeneRanges(data_frame(gene_id = "ATAC_peak_261927"), atac_list$gene_metadata, 1e5)
peak1_rasqual_pvalues = purrr::map_df(qtlResults()$atac_rasqual, ~tabixFetchGenes(peak1_region, .)[[1]], .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344"))) %>%
  dplyr::arrange(p_nominal) %>%
  dplyr::group_by(condition_name) %>%
  dplyr::mutate(R2 = calculateR2FromLead(snp_id, vcf_file$genotypes)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(chr = as.character(chr))

peak1_plot = ggplot(peak1_rasqual_pvalues, aes(x = pos, y = -log(p_nominal, 10), color = R2)) + 
  geom_point() + 
  facet_grid(condition_name~., scales = "free_y") +
  theme_light()
ggsave("figures/main_figures/PTK2B_peak1_rasqual_manhattan.pdf", plot = peak1_plot, width = 6, height = 7)


#Import Alzheimer's GWAS summary stats
igap_gwas_results = readRDS("annotations/IGAP_stage_1_2_combined_stats.rds")
gwas_pvalues = dplyr::filter(igap_gwas_results, chr == 8, pos > min(joint_pvalues$pos), pos < max(joint_pvalues$pos)) %>%
  dplyr::transmute(chr = as.integer(chr), pos, p_nominal = igap_pvalue, beta = Beta, se = SE, condition = "AD GWAS")
gwas_pvalues_df = dplyr::select(naive_pvalues_linear, chr, pos, snp_id, MAF) %>% 
  dplyr::left_join(gwas_pvalues, ., by = c("chr", "pos"))

#Test for colocalisation
a = coloc.abf(dataset1 = list(pvalues = naive_pvalues_linear$p_nominal, N = 84, MAF = naive_pvalues_linear$MAF, 
                              type = "quant", beta = naive_pvalues_linear$beta, snp = naive_pvalues_linear$snp_id), 
              dataset2 = list(pvalues = gwas_pvalues_df$p_nominal, MAF = gwas_pvalues_df$MAF, beta = gwas_pvalues_df$beta, 
                              snp = gwas_pvalues_df$snp_id, N = 74046, type = "cc", s = 0.35))  



pvalues_df = dplyr::select(naive_pvalues_linear, chr, pos, snp_id, MAF) %>% 
  dplyr::left_join(pvalues, ., by = c("chr", "pos")) %>%
  dplyr::filter(!is.na(snp_id))

a = coloc.abf(dataset1 = list(pvalues = naive_pvalues_linear$p_nominal, N = 84, MAF = naive_pvalues_linear$MAF, 
                              type = "quant", beta = naive_pvalues_linear$beta, snp = naive_pvalues_linear$snp_id), 
              dataset2 = list(pvalues = pvalues_df$p_nominal, MAF = pvalues_df$MAF, beta = pvalues_df$OR, 
                              snp = pvalues_df$snp_id, N = 74046, type = "cc", s = 0.35))  


#Plot all pvalues
all_pvalues = rbind(dplyr::select(joint_pvalues, pos, p_nominal, condition), gwas_pvalues)
all_pvalues$condition = factor(all_pvalues$condition, levels = c("AD GWAS", "naive", "IFNg_SL1344"))
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/eQTL_vs_GWAS_manhattan_plot.pdf", plot = eQTL_manhattan_plot, width = 7, height = 8)

#Make effect size plots for the ATAC QTLs
#Make effect size plots for these two qtls
naive_QTL = plotEQTL("ATAC_peak_261927", "rs28834970", atac_list$cqn, vcf_file$genotypes, 
                      atac_list$sample_metadata, atac_list$gene_metadata)
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/PTK2B_enhancer1_QTL.pdf", naive_QTL, width = 7, height = 7)

ifng_sl1344_QTL = plotEQTL("ATAC_peak_261942", "rs1429938", atac_list$cqn, vcf_file$genotypes, 
                           atac_list$sample_metadata, atac_list$gene_metadata)
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/PTK2B_enhancer2_eQTL.pdf", ifng_sl1344_QTL, width = 7, height = 7)



#Test for colocalisation
a = testColoc(naive_pvalues, naive_atac_pvalues, n1 = 69, n2 = 42)
b = testColoc(naive_pvalues, IFNg_SL1344_atac_pvalues, n1 = 69, n2 = 42)
c = testColoc(IFNg_SL1344_pvalues, IFNg_SL1344_atac_pvalues, n1 = 69, n2 = 31)
c = testColoc(IFNg_SL1344_pvalues, naive_pvalues, n1 = 69, n2 = 69)



#Perform finemapping with ATAC data

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Overlap
naive_credible_set = constructCredibleSet(naive_pvalues, threshold = 0.99)
naive_cs_annotated = annotateCredibleSet(naive_credible_set, atac_list$gene_metadata, atac_tabix_list$naive)

ifng_sl1344_cs = constructCredibleSet(IFNg_SL1344_pvalues, threshold = 0.99)
ifng_sl1344_cs_annotated = annotateCredibleSet(ifng_sl1344_cs, atac_list$gene_metadata, atac_tabix_list$IFNg_SL1344)



#Check if the causal SNPs disrupt any known TF motifs
#Import motif matches
motif_metadata = readRDS("results/ATAC/cisBP/cisBP_motif_metadata.rds") %>%
  dplyr::transmute(motif_id = Motif_ID, tf_name = TF_Name, tf_count = TF_count)
motif_disruptions = importMotifDisruptions("results/ATAC/motif_analysis/motif_disruption.txt") %>%
  dplyr::left_join(motif_metadata, by = "motif_id")

#Filter by SNP ID
motif_hits = dplyr::filter(motif_disruptions, snp_id %in% causal_variants) %>%
  dplyr::filter(max_rel_score > 0.8) %>% dplyr::arrange(-abs(rel_diff))

#PTK2B and CLU hits are independent from each other
calculateR2FromLead(c("rs7982", "rs28834970"), vcf_file$genotypes)


#Find associated peaks for a single specific SNP
atac_snp_results = purrr::map(qtlResults()$atac_rasqual, ~tabixFetchSNPsQuick(c("rs11375466"),.,vcf_file$snpspos)) %>%
  ldply(.id = "condition_name")
dplyr::filter(atac_snp_results, p_nominal < 1e-4) %>% 
  dplyr::select(gene_id, condition_name, p_nominal) %>% 
  dplyr::left_join(dplyr::select(atac_data$gene_metadata, gene_id, start, end))
