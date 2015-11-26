library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("plyr")
library("qvalue")

#Import expression data
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")

#Import genotypes
vcf_file = readRDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")

#Import eQTL mapping hits
fastqtl_callset = readRDS("results/SL1344/fastqtl/output/fastqtl_call_set.rds")
interaction_snps = unique(fastqtl_callset$interaction$snp_id)

#Construct genepos and SNP-pos matrices
atac_snpspos = dplyr::filter(vcf_file$snpspos, snpid %in% interaction_snps) %>% 
  as.data.frame()
atac_genepos = dplyr::rename(atac_list$gene_metadata, geneid = gene_id) %>% 
  dplyr::select(geneid, chr, left, right)
IFNg_atac_hits = runMatrixEQTL(exp_data = atac_list$exprs_cqn_list$IFNg, 
                    geno_data = extractSubset(atac_list$design_list$IFNg, vcf_file$genotypes, "genotype_id", "donor"),
                    snpspos = atac_snpspos, 
                    genepos = atac_genepos, 
                    cisDist = 2e05, pvOutputThreshold = 1)

naive_atac_hits = runMatrixEQTL(exp_data = atac_list$exprs_cqn_list$naive, 
                               geno_data = extractSubset(atac_list$design_list$naive, vcf_file$genotypes, "genotype_id", "donor"),
                               snpspos = atac_snpspos, 
                               genepos = atac_genepos, 
                               cisDist = 2e05, pvOutputThreshold = 1)

#Test for interaction betweeb IFNg and naive conditions
simple_model <- function(model_data){
  model = lm(expression ~ genotype + condition_name + sex + ng_ul_mean + diff_days + PEER_factor_1 + PEER_factor_2 + PEER_factor_3 + PEER_factor_4, model_data)
  return(model)
}
interaction_model <- function(model_data){
  model = lm(expression~genotype + condition_name + condition_name:genotype + sex + ng_ul_mean + diff_days + PEER_factor_1 + PEER_factor_2 + PEER_factor_3 + PEER_factor_4, model_data)
  return(model)
}

ifng_eqtl_data_list = eqtl_data_list
ifng_eqtl_data_list$sample_metadata = dplyr::filter(ifng_eqtl_data_list$sample_metadata, condition %in% c("A","B"))
ifng_eqtl_data_list$exprs_cqn = ifng_eqtl_data_list$exprs_cqn[,ifng_eqtl_data_list$sample_metadata$sample_id]
ifng_eqtl_data_list$covariates = ifng_eqtl_data_list$covariates[ifng_eqtl_data_list$sample_metadata$sample_id,]

interaction_pvalues = testMultipleInteractions(fastqtl_callset$interaction, ifng_eqtl_data_list, vcf_file, simple_model, interaction_model, return = "all")
pvals = lapply(interaction_pvalues, function(x){p = x$anova[[6]][2]})

interaction_df = plyr::ldply(pvals, .id = "id") %>% 
  tidyr::separate(id, into = c("gene_id", "snp_id"), sep = ":") %>%
  dplyr::rename(pvalue = V1) %>% 
  dplyr::left_join(gene_id_name_map, by = "gene_id") %>%
  tbl_df() %>%
  dplyr::select(gene_id, gene_name, snp_id, pvalue) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::mutate(qvalue = qvalue(pvalue)$qvalues)

ifng_interactions = dplyr::filter(interaction_df, qvalue < 0.1)

#Perform the same analysis for RNA-Seq
#Construct genepos and SNP-pos matrices
rna_snpspos = dplyr::filter(vcf_file$snpspos, snpid %in% interaction_snps) %>% 
  as.data.frame()
rna_genepos = dplyr::rename(eqtl_data_list$gene_metadata, geneid = gene_id, chr = chromosome_name, left = start_position, right = end_position) %>% 
  dplyr::select(geneid, chr, left, right) %>% as.data.frame()
IFNg_rna_hits = runMatrixEQTL(exp_data = eqtl_data_list$exprs_cqn_list$IFNg, 
                              geno_data = extractSubset(dplyr::filter(eqtl_data_list$sample_metadata, condition_name == "IFNg"), 
                                                        vcf_file$genotypes, "genotype_id", "donor"),
                              snpspos = rna_snpspos, 
                              genepos = rna_genepos, 
                              cisDist = 5e05, pvOutputThreshold = 1)

naive_rna_hits = runMatrixEQTL(exp_data = eqtl_data_list$exprs_cqn_list$naive, 
                               geno_data = extractSubset(dplyr::filter(eqtl_data_list$sample_metadata, condition_name == "naive"), 
                                                         vcf_file$genotypes, "genotype_id", "donor"),
                               snpspos = rna_snpspos, 
                               genepos = rna_genepos, 
                               cisDist = 5e05, pvOutputThreshold = 1)

IFNg_rna_betas = dplyr::rename(IFNg_rna_hits$cis$eqtls, snp_id = snps, gene_id = gene, ifng_beta = beta, ifng_pvalue = pvalue) %>% 
  dplyr::semi_join(ifng_interactions, by = c("snp_id", "gene_id"))
naive_rna_betas = dplyr::rename(naive_rna_hits$cis$eqtls, snp_id = snps, gene_id = gene, naive_beta = beta, naive_pvalue = pvalue) %>% 
  dplyr::semi_join(ifng_interactions, by = c("snp_id", "gene_id"))
rna_joint_effects = dplyr::left_join(IFNg_rna_betas, naive_rna_betas, by = c("snp_id", "gene_id"))
rna_joint_effects_2 = rna_joint_effects[which(abs(abs(rna_joint_effects$naive_beta) - abs(rna_joint_effects$ifng_beta)) > 0.2),]

ifng_interaction_plot = ggplot(rna_joint_effects_2, aes(x = naive_beta, y = ifng_beta)) + geom_point() + 
  xlab("Effect size in naive condition") + 
  ylab("Effect size in IFNg condition")
ggsave("results/SL1344/eQTLs/RNA_effect_size_plot.pdf", ifng_interaction_plot, width = 5, height = 5)


#Join p-values and effect sizes from both conditions
IFNg_lead_peaks = dplyr::group_by(IFNg_atac_hits$cis$eqtls, snps) %>% 
  arrange(pvalue) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::rename(gene_id = gene, snp_id = snps) %>% 
  dplyr::select(gene_id, snp_id) %>% ungroup() %>% 
  dplyr::semi_join(rna_join_effects_2, by = "snp_id")
naive_lead_peaks = dplyr::group_by(naive_atac_hits$cis$eqtls, snps) %>% 
  arrange(pvalue) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::rename(gene_id = gene, snp_id = snps) %>% 
  dplyr::select(gene_id, snp_id) %>% ungroup() %>% 
  dplyr::semi_join(rna_join_effects_2, by = "snp_id")
joint_peaks = rbind(naive_lead_peaks, IFNg_lead_peaks) %>% unique()

naive_betas = dplyr::rename(naive_atac_hits$cis$eqtls, gene_id = gene, snp_id = snps) %>%
  dplyr::semi_join(joint_peaks, by = c("gene_id", "snp_id")) %>% 
  dplyr::rename(naive_beta = beta, naive_pvalue = pvalue)

IFNg_betas = dplyr::rename(IFNg_atac_hits$cis$eqtls, gene_id = gene, snp_id = snps) %>%
  dplyr::semi_join(joint_peaks, by = c("gene_id", "snp_id")) %>% 
  dplyr::rename(ifng_beta = beta, ifng_pvalue = pvalue)

atac_joint_effects = dplyr::left_join(naive_betas, IFNg_betas, by = c("snp_id", "gene_id"))

#Look at their effect sizes in two conditions
atac_plot = ggplot(atac_joint_effects, aes(x = naive_beta, y = ifng_beta)) + geom_point() + 
  xlab("Effect size in naive condition") +
  ylab("Effect size in IFNg condition")
ggsave("results/SL1344/eQTLs/atac_effect_size_plot.pdf", atac_plot, width = 5, height = 5)

