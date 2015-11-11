library("devtools")
library("qvalue")
library("dplyr")
library("ggplot2")
library("SNPRelate")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Import data
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Import genotypes
SNPRelate::snpgdsVCF2GDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz", 
                         "results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.gds", method = "copy.num.of.ref")
vcf_file = gdsToMatrix("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.gds")
saveRDS(vcf_file, "results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")
vcf_file = readRDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")

#Import permutation p-values
naive_qtls = importFastQTLTable("results/SL1344/fastqtl/output/naive_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
IFNg_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/SL1344_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
IFNg_SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_SL1344_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)

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

#Merge conditions together
naive_hits = dplyr::filter(naive_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, qvalue)
IFNg_hits = dplyr::filter(IFNg_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, qvalue)
SL1344_hits = dplyr::filter(SL1344_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, qvalue)
IFNg_SL1344_hits = dplyr::filter(IFNg_SL1344_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, qvalue)

#Merge all of the hits together
joint_df = dplyr::full_join(naive_hits, IFNg_hits, by = c("gene_id","snp_id")) %>%
  dplyr::full_join(SL1344_hits, c("gene_id","snp_id")) %>%
  dplyr::full_join(IFNg_SL1344_hits,c("gene_id","snp_id"))

#Test for interaction
interaction_pvalues = testMultipleInteractions(joint_df, eqtl_data_list, vcf_file)
saveRDS(interaction_pvalues, "results/SL1344/interaction_pvalues.rds")

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
plotEQTL("ENSG00000091490", "rs13141111", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)
plotEQTL("ENSG00000091490", "rs6831024", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)

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


