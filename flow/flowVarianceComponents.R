library("dplyr")
library("lme4")
library("devtools")
load_all("../seqUtils/")
library("GWASTools")
library("ggplot2")
library("MatrixEQTL")

#### Functions ####
mapFlowQTLs <- function(intensity_df, genepos, genotype_list, genotype_donor_map, cisDist = 5e5, covariates = NULL, permute = FALSE){
  #Helper function to map QTLs for flow data
  genotype_matrix = extractSubset(genotype_donor_map, genotype_list$genotypes, 
                                  old_column_names = "genotype_id", new_column_names = "donor")
  print(head(genotype_matrix))
  qtl_results = runMatrixEQTL(intensity_df, genotype_matrix, as.data.frame(genotype_list$snpspos), 
                              genepos, cisDist, pvOutputThreshold = 1, covariates = covariates, permute)
  qtl_df = dplyr::left_join(qtl_results$cis$eqtls, genotype_list$snpspos, by = c("snps" = "snpid")) %>%
    dplyr::mutate(log10_pvalue = -log(pvalue,10)) %>%
    dplyr::rename(snp_id = snps) %>%
    dplyr::mutate(expected = -log(c(1:length(pvalue))/length(pvalue),10))
  return(qtl_df)
}


#Import data
flow_purity = readRDS("macrophage-gxe-study/data/covariates/flow_cytometry_purity.rds")
line_medatada = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")
vcf_file = readRDS("genotypes/SL1344/array_genotypes.59_samples.imputed.vcfToMatrix.rds")
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Extract gene coords
gene_coords = dplyr::filter(eqtl_data_list$genepos, geneid %in% 
                              c("ENSG00000170458","ENSG00000203747","ENSG00000162747","ENSG00000260314")) %>%
  dplyr::left_join(gene_id_name_map, by = c("geneid" = "gene_id"))

#Map channel to marker
channel_marker_map = data_frame(channel = c("APC.A","PE.A","Pacific.Blue.A"), marker = c("CD206","CD16","CD14"))
unique_lines = dplyr::select(line_medatada, line_id, donor, genotype_id) %>% unique()
flow_data = dplyr::left_join(flow_purity, channel_marker_map, by = "channel") %>%
  dplyr::left_join(unique_lines, by = "donor") %>%
  dplyr::mutate(donor = ifelse(donor == "fpdj", "nibo",donor)) %>% #fpdj and nibo are the same donors
  dplyr::mutate(intensity = mean2-mean1) %>%
  dplyr::select(line_id, genotype_id, donor, flow_date, marker, purity, intensity)

#Import genotypes
#Cd14
#SNPRelate::snpgdsVCF2GDS("CD14_cis_region.vcf", "CD14_cis_region1.gds",method = "copy.num.of.ref")
cd14_cis_region = gdsToMatrix("CD14_cis_region1.gds")

#CD16 (FCGR3B)
#bcftools view -r 1:161041759-162131963 -S macrophage-gxe-study/data/sample_lists/flow_cytometry_gt_list.txt genotypes/GRCh38/imputed_20151005/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr1.GRCh38.sorted.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - | bcftools norm -m+both - | bcftools view -m2 -M2 - > CD16_cis.vcf
#SNPRelate::snpgdsVCF2GDS("CD16_cis.vcf", "CD16_cis.gds",method = "copy.num.of.ref")
CD16_cis_region = gdsToMatrix("CD16_cis.gds")

#### Estimate variance explained by donor ####
#Analyse CD14 only
cd14_geno_mat = as.data.frame(t(cd14_cis_region$genotypes)) %>% 
  dplyr::mutate(genotype_id = colnames(cd14_cis_region$genotypes))
cd14_data = dplyr::filter(flow_data, marker == "CD14")
cd14_model_data = dplyr::left_join(cd14_data, cd14_geno_mat, by = "genotype_id")
cd14_model = lmer(intensity ~ (1|flow_date) + (1|donor), cd14_model_data)
cd14_variance = seqUtils::varianceExplained(cd14_model)

#Analyse CD16 only
cd16_geno_mat = as.data.frame(t(CD16_cis_region$genotypes)) %>% 
  dplyr::mutate(genotype_id = colnames(CD16_cis_region$genotypes))
cd16_data = dplyr::filter(flow_data, marker == "CD16")
cd16_model_data = dplyr::left_join(cd16_data, cd16_geno_mat, by = "genotype_id")
cd16_model = lmer(intensity ~ (1|flow_date) + (1|donor), cd16_model_data, REML = TRUE)
cd16_variance = seqUtils::varianceExplained(cd16_model)

#Analyse CD206 only
cd206_data = dplyr::filter(flow_data, marker == "CD206")
cd206_model = lmer(intensity ~ (1|flow_date) + (1|donor), cd206_data)
cd206_variance = seqUtils::varianceExplained(cd206_model)

variance_explained = dplyr::mutate(rbind(cd14_variance, cd16_variance, cd206_variance), marker = c("CD14","CD16","CD206")) %>% 
  dplyr::select(marker, everything(), -type)
write.table(variance_explained, "results/flow/VarComp/variance_explained.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#### QTL mapping ####
##### Make a single intensity data set #####
cd14_mean = dplyr::group_by(cd14_data, donor) %>% summarise(CD14 = mean(intensity))
cd16_mean = dplyr::group_by(cd16_data, donor) %>% summarise(CD16 = mean(intensity))
cd206_mean = dplyr::group_by(cd206_data, donor) %>% summarise(CD206 = mean(intensity))

#Prepare intensity data
intensity_df = dplyr::left_join(cd14_mean, cd16_mean, by = "donor") %>% 
  dplyr::left_join(cd206_mean, by = "donor") %>% 
  as.data.frame()
rownames(intensity_df) = intensity_df$donor
intensity_df = t(intensity_df[-2,-1])

#Prepare gene positions
genepos = dplyr::select(gene_coords, gene_name, chr, left, right) %>%
  dplyr::rename(geneid = gene_name)
genepos$geneid = c("CD14","CD206", "CD16","CD16")

#Prepare genotype data
genotype_donor_map = dplyr::select(flow_data, genotype_id, donor) %>% unique() %>% 
  dplyr::filter(genotype_id != "HPSI0813i-fpdj_3")

#### CD14 QTLs ####
cd14_qtl_df = mapFlowQTLs(intensity_df, genepos[1:3,], cd14_cis_region, genotype_donor_map, cisDist = 5e5)
cd14_manhattan = ggplot(cd14_qtl_df, aes(x = pos, y = log10_pvalue)) + geom_point()
ggsave("results/flow/VarComp/CD14_manhattan_plot.pdf", plot = cd14_manhattan, width = 8, height = 6)

#Make a QQ plot
ggplot(cd14_qtl_df, aes(x = expected, y = log10_pvalue)) + 
  geom_point() + 
  stat_abline(slope = 1, intercept = 0, color = "red")

#Make a QTL plot
cd14_qtl_plot = ggplot(cd14_model_data, aes(x = factor(rs2569177), y = intensity)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = .1))
ggsave("results/flow/VarComp/CD14_lead_QTL.pdf", plot = cd14_qtl_plot, width = 6, height = 6)

#Perform VarComp analysis with the lead SNP
cd14_model_qtl = lmer(intensity ~ (1|flow_date) + (1|donor) + (1|rs2569177), cd14_model_data)
cd14_variance_qtl = seqUtils::varianceExplained(cd14_model_qtl)
write.table(cd14_variance_qtl, "results/flow/VarComp/CD14_variance_explained.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#### CD16 QTLS ####
cd16_qtl_df = mapFlowQTLs(intensity_df, genepos[1:3,], CD16_cis_region, genotype_donor_map, cisDist = 2e5)
cd16_manhattan = ggplot(cd16_qtl_df, aes(x = pos, y = log10_pvalue)) + geom_point()
ggsave("results/flow/VarComp/CD16_manhattan_plot.pdf", plot = cd16_manhattan, width = 8, height = 6)

#Make a QQ plot
ggplot(cd16_qtl_df, aes(x = expected, y = log10_pvalue)) + 
  geom_point() + 
  stat_abline(slope = 1, intercept = 0, color = "red")

#Perform permutations
perm_list = lapply(as.list(c(1:100)), function(x){mapFlowQTLs(intensity_df, genepos[1:3,], CD16_cis_region, genotype_donor_map, cisDist = 2e5, permute = TRUE)})
perm_max_statistic = lapply(perm_list, function(x){max(abs(x$statistic))})


#Fit VarComp model again with four independent associations
cd16_model_qtl = lmer(intensity ~ (1|flow_date) + (1|donor) + (1|rs10917809) + (1|rs4657019) + (1|rs2333845) + (1|rs4571943), cd16_model_data, REML = TRUE)
cd16_variance_qtl = seqUtils::varianceExplained(cd16_model_qtl)
anova(cd16_model_qtl, cd16_model)
write.table(cd16_variance_qtl, "results/flow/VarComp/CD16_variance_explained.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Test these four SNPs in a linear model
cd16_lm = lm(intensity ~ 1,cd16_model_data)
cd16_snps = lm(intensity ~ rs10917809 + rs2333845 + rs4657019 + rs4571943, cd16_model_data)
anova(cd16_lm, cd16_snps)

#Make sure that the SNPs are independent
cd16_qtl_snps = c("rs10917809","rs2333845","rs4657019")
cd16_qtl_snp_rsquared = cor(t(CD16_cis_region$genotypes[cd16_qtl_snps,]))^2
write.table(cd16_qtl_snp_rsquared, "results/flow/VarComp/CD16_qtl_snp_rsquared.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Make boxplots of these SNPs
cd16_rs4657019_plot = ggplot(cd16_model_data, aes(x = factor(rs4657019), y = intensity)) + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(width = .1))
ggsave("results/flow/VarComp/CD16_rs4657019_plot.pdf", plot = cd16_rs4657019_plot, width = 6, height = 6)
cd16_rs10917809_plot = ggplot(cd16_model_data, aes(x = factor(rs10917809), y = intensity)) +  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(width = .1))
ggsave("results/flow/VarComp/CD16_rs10917809_plot.pdf", plot = cd16_rs10917809_plot, width = 6, height = 6)
cd16_rs2333845_plot = ggplot(cd16_model_data, aes(x = factor(rs2333845), y = intensity)) +  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(width = .1))
ggsave("results/flow/VarComp/CD16_rs2333845_plot.pdf", plot = cd16_rs2333845_plot, width = 6, height = 6)
cd16_rs4571943_plot = ggplot(cd16_model_data, aes(x = factor(rs4571943), y = intensity)) +  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(width = .1))
ggsave("results/flow/VarComp/CD16_rs4571943_plot.pdf", plot = cd16_rs4571943_plot, width = 6, height = 6)

#Look at CD16 gene expression
exp_mat = t(eqtl_data_list$exprs_cqn_list$naive[c("ENSG00000203747","ENSG00000162747"),]) 
donor_genotype_map = dplyr::filter(eqtl_data_list$sample_metadata, condition == "A") %>% dplyr::select(donor, genotype_id)
cd16_exp = dplyr::mutate(as.data.frame(exp_mat), donor = rownames(exp_mat)) %>% 
  dplyr::rename(FCGR3A = ENSG00000203747, FCGR3B = ENSG00000162747) %>% 
  dplyr::left_join(donor_genotype_map, by = "donor") %>%
  dplyr::left_join(cd16_geno_mat, by = "genotype_id") %>%
  dplyr::filter(!is.na(rs4657019))

fcgr3a_rs4657019_plot = ggplot(cd16_exp, aes(x = factor(rs4657019), y = FCGR3A)) +  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(width = .1))
ggsave("results/flow/VarComp/FCGR3A_rs4657019_plot.pdf", plot = fcgr3a_rs4657019_plot, width = 6, height = 6)
fcgr3a_rs2333845_plot = ggplot(cd16_exp, aes(x = factor(rs2333845), y = FCGR3A)) +  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(width = .1))
ggsave("results/flow/VarComp/FCGR3A_rs2333845_plot.pdf", plot = fcgr3a_rs2333845_plot, width = 6, height = 6)
fcgr3b_rs10917809_plot = ggplot(cd16_exp, aes(x = factor(rs10917809), y = FCGR3B)) +  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(width = .1))
ggsave("results/flow/VarComp/FCGR3B_rs10917809_plot.pdf", plot = fcgr3b_rs10917809_plot, width = 6, height = 6)
