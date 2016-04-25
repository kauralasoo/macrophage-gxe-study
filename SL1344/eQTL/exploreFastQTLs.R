library("devtools")
library("qvalue")
library("dplyr")
library("ggplot2")
library("SNPRelate")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import genotypes
#SNPRelate::snpgdsVCF2GDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz", 
#                         "results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.gds", method = "copy.num.of.ref")
#vcf_file = gdsToMatrix("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.gds")
#saveRDS(vcf_file, "results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")
#vcf_file = readRDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import permutation p-values
naive_qtls = importFastQTLTable("results/SL1344/fastqtl/output/naive_500kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
IFNg_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_500kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/SL1344_500kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
IFNg_SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_SL1344_500kb_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)

#Combine all hits to a list
qtl_list = list(naive = naive_qtls, IFNg = IFNg_qtls, SL1344 = SL1344_qtls, IFNg_SL1344 = IFNg_SL1344_qtls)
pi1_matrix = calculatePairwisePi1(qtl_list)
write.table(pi1_matrix, "results/SL1344/eQTLs/pi1_matrix.txt", quote = FALSE, sep = "\t")
pi1_matrix_tidy = calculatePairwisePi1(qtl_list, tidy = TRUE)

#Make a QQ plot
ggplot(naive_qtls, aes(x = p_expected, y = p_beta_log10)) +
  geom_point() +
  stat_abline(slope = 1, intercept = 0, color = "red")

#Make some example plots
cd14_plot = plotEQTL("ENSG00000170458", "rs778583", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
                     eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)
ctsc_plot = plotEQTL("ENSG00000109861", "rs11019479", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
                     eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)
ccs_plot = plotEQTL("ENSG00000173992", "rs636128", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
                     eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)

#Merge conditions together
naive_hits = dplyr::filter(naive_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, p_nominal, p_beta, qvalue)
IFNg_hits = dplyr::filter(IFNg_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, p_nominal, p_beta, qvalue)
SL1344_hits = dplyr::filter(SL1344_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, p_nominal, p_beta, qvalue)
IFNg_SL1344_hits = dplyr::filter(IFNg_SL1344_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, p_nominal, p_beta, qvalue)

#Merge all of the hits together
joint_df = dplyr::full_join(naive_hits, IFNg_hits, by = c("gene_id","snp_id")) %>%
  dplyr::full_join(SL1344_hits, c("gene_id","snp_id")) %>%
  dplyr::full_join(IFNg_SL1344_hits,c("gene_id","snp_id"))

#Test for interaction
simple_model <- function(model_data){
  model = lm(expression ~ genotype + condition_name + sex + ng_ul_mean + diff_days + PEER_factor_1 + PEER_factor_2 + PEER_factor_3 + PEER_factor_4, model_data)
  return(model)
}
interaction_model <- function(model_data){
  model = lm(expression~genotype + condition_name + condition_name:genotype + sex + ng_ul_mean + diff_days + PEER_factor_1 + PEER_factor_2 + PEER_factor_3 + PEER_factor_4, model_data)
  return(model)
}

interaction_pvalues = testMultipleInteractions(joint_df, eqtl_data_list, vcf_file, simple_model, interaction_model)
saveRDS(interaction_pvalues, "results/SL1344/interaction_pvalues.rds")
interaction_pvalues = readRDS("results/SL1344/interaction_pvalues.rds")

#Reformat the interaction data frame
interaction_df = plyr::ldply(interaction_pvalues, .id = "id") %>% 
  tidyr::separate(id, into = c("gene_id", "snp_id"), sep = ":") %>%
  dplyr::rename(pvalue = V2) %>% 
  dplyr::left_join(gene_id_name_map, by = "gene_id") %>%
  tbl_df() %>%
  dplyr::select(gene_id, gene_name, snp_id, pvalue) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::mutate(qvalue = qvalue(pvalue)$qvalues)

#Keep the stronges interaction per gene
interaction_filtered = dplyr::filter(interaction_df, qvalue < 0.1) %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(gene_id, qvalue) %>% 
  dplyr::top_n(1) %>% 
  dplyr::ungroup() %>% dplyr::arrange(qvalue)

#Make a list of all QTLs
qtl_list = list(naive = naive_hits, IFNg = IFNg_hits, SL1344 = SL1344_hits, 
                IFNg_SL1344 = IFNg_SL1344_hits, interaction = interaction_filtered)
saveRDS(qtl_list, "results/SL1344/eQTLs/fastqtl_call_set.rds")



plotEQTL("ENSG00000196735", "rs79617990", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)

interaction_plots = makeMultiplePlots(interaction_df[1:200,],eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
                                      eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)

savePlots(interaction_plots, "results/SL1344/eQTLs/interaction_plots/", width = 7, height = 7)

multi_snp_genes = group_by(joint_df[,1:2], gene_id) %>% 
  dplyr::summarise(snp_count = length(gene_id)) %>% 
  dplyr::filter(snp_count == 2)


snps_per_gene = dplyr::semi_join(joint_df[,1:2], multi_snp_genes, by = "gene_id") %>%
  dplyr::select(gene_id, snp_id) %>%
  group_by(gene_id) %>% 
  dplyr::summarize(snps = paste(snp_id, collapse = "€")) %>% 
  tidyr::separate(snps, c("snp1","snp2"), sep ="€")

#Calculate LD for each pair
snps = union(snps_per_gene$snp1, snps_per_gene$snp2)
selected_snps = vcf_file$genotypes[snps, ]
ld_mat = dplyr::group_by(snps_per_gene, gene_id) %>% 
  dplyr::mutate(ld = snpgdsLDpair(selected_snps[snp1,],selected_snps[snp2,], method = "r")[1])
dplyr::filter(ld_mat, ld*ld < 0.2) %>% View()

ld_mat2 = dplyr::group_by(snps_per_gene, gene_id) %>% 
  dplyr::mutate(ld = snpgdsLDpair(selected_snps[snp1,],selected_snps[snp2,], method = "dprime")[1])
dplyr::filter(ld_mat2, ld < 0.8)

plot(ld_mat$ld*ld_mat$ld, ld_mat2$ld*ld_mat2$ld)


#Regress on lead SNP to find alternatives
#Example SEL1L3 gene
plotEQTL("ENSG00000091490", "rs10006078", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
plotEQTL("ENSG00000091490", "rs6831024", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)

#Fetch p-values from rasqual output
gene_df = data_frame(gene_id = "ENSG00000091490")
gene_ranges = constructGeneRanges(gene_df, combined_expression_data$gene_metadata, cis_window = 5e5)
tabix_data = tabixFetchGenes(gene_ranges, "databases/SL1344/IFNg_SL1344_500kb.sorted.txt.gz")
cond1 = tabix_data[[1]] %>% dplyr::filter(p_nominal < 0.01)
plot(a$pos, a$chisq)

gene_df = data_frame(gene_id = "ENSG00000091490")
gene_ranges = constructGeneRanges(gene_df, combined_expression_data$gene_metadata, cis_window = 5e5)
tabix_data = tabixFetchGenes(gene_ranges, "databases/SL1344/IFNg_500kb.sorted.txt.gz")
cond2 = tabix_data[[1]] %>% dplyr::filter(p_nominal < 0.01)
plot(a$pos, a$chisq)
fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }
cond2$z = fisher.r2z(sqrt(cond2$chisq/69)*sign(cond2$beta))

#Export effect sizes
write.table(dplyr::select(cond1, snp_id, beta), "SEL1L3/ifng_sl1344_log2FC.z", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
write.table(dplyr::select(cond1, snp_id, beta), "SEL1L3/ifng_log2FC.z", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

snp_correlation_matrix = cor(t(vcf_file$genotypes[cond1$snp_id,]))
write.table(snp_correlation_matrix, "SEL1L3/ifng_log2FC.ld", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

#Set up covariates
covariates = eqtl_data_list$covariates_list$IFNg_SL1344[1:7,]
cov_mat = dplyr::mutate(as.data.frame(t(covariates)), donor = colnames(covariates))

expression = eqtl_data_list$exprs_cqn_list$IFNg_SL1344["ENSG00000091490",]
exp_mat = data_frame(expression = expression, donor = names(expression))

#Sort out genotypes
genotypes = vcf_file$genotypes[c("rs6831024","rs13141111"),] %>% t() %>% as.data.frame() %>%
  dplyr::mutate(genotype_id = colnames(vcf_file$genotypes))
geno_donor_map = eqtl_data_list$sample_metadata %>% dplyr::select(donor, genotype_id) %>% unique()

df = dplyr::left_join(exp_mat, cov_mat, by = "donor") %>%
  dplyr::left_join(geno_donor_map, by = "donor") %>%
  dplyr::left_join(genotypes, by = "genotype_id")

#Calculate residuals
m1 = lm(expression ~ sex + ng_ul_mean + diff_days + PEER_factor_1 + PEER_factor_2 + PEER_factor_3 + PEER_factor_4 + rs13141111, data = df)
df$res = m1$residuals

ggplot(df, aes(x = factor(rs13141111), y = res)) + geom_boxplot() + geom_point()
ggplot(df, aes(x = factor(rs6831024), y = res)) + geom_boxplot() + geom_point()


#Plot IL2RA expression
plotEQTL("ENSG00000134460", "rs7911500", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)
plotEQTL("ENSG00000134460", "rs12722605", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)


