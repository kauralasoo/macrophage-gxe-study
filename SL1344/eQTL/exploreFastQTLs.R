library("devtools")
library("qvalue")
library("dplyr")
load_all("../seqUtils/")

#Helper functions
enrichFastQTLPvalues <- function(fastqtl_pvalues, gene_id_name_map){
  res = tbl_df(fastqtl_pvalues) %>%
    dplyr::arrange(p_beta) %>%
    dplyr::filter(!is.na(p_beta)) %>%
    dplyr::mutate(p_beta_log10 =-log(p_beta,10)) %>% 
    dplyr::mutate(p_expected = -log(c(1:length(p_beta))/length(p_beta),10)) %>% 
    dplyr::left_join(gene_id_name_map, by = "gene_id") %>%
    dplyr::select(gene_id, gene_name, everything()) %>%
    dplyr::mutate(qvalue = qvalue::qvalue(p_beta)$qvalues)
  return(res)
}

#Import data
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Import genotypes
#SNPRelate::snpgdsVCF2GDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.vcf.gz", 
                         "results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.gds", method = "copy.num.of.ref")
#vcf_file = gdsToMatrix("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.gds")
#saveRDS(vcf_file, "results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.rds")
vcf_file = readRDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.rds")

#Import permutation p-values
naive_qtls = importFastQTLTable("results/SL1344/fastqtl/output/naive_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
IFNg_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/SL1344_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
IFNg_SL1344_qtls = importFastQTLTable("results/SL1344/fastqtl/output/IFNg_SL1344_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)

#Combine all hits to a list
qtl_list = list(naive = naive_qtls, IFNg = IFNg_qtls, SL1344 = SL1344_qtls, IFNg_SL1344 = IFNg_SL1344_qtls)
pi1_matrix = calculatePairwisePi1(qtl_list)
pi1_matrix_tidy = calculatePairwisePi1(qtl_list, tidy = TRUE)

#Make a QQ plot
ggplot(IFNg_SL1344_qtls, aes(x = p_expected, y = p_beta_log10)) +
  geom_point() +
  stat_abline(slope = 1, intercept = 0, color = "red")

#Make some example plots
cd14_plot = plotEQTL("ENSG00000170458", "rs778583", eqtl_dataset$exprs_cqn, vcf_file$genotypes, 
                     eqtl_dataset$sample_metadata, eqtl_dataset$gene_metadata)


#Merge conditions together
naive_hits = dplyr::filter(naive_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, qvalue)
IFNg_hits = dplyr::filter(IFNg_qtls, qvalue < 0.1) %>% dplyr::select(gene_id, snp_id, qvalue)
joint_df = dplyr::full_join(naive_hits, IFNg_hits, by = c("gene_id","snp_id"))

multi_snp_genes = group_by(joint_df, gene_id) %>% 
  dplyr::summarise(snp_count = length(gene_id)) %>% 
  dplyr::filter(snp_count == 2)


snps_per_gene = dplyr::semi_join(joint_df, multi_snp_genes, by = "gene_id") %>%
  dplyr::select(gene_id, snp_id) %>%
  group_by(gene_id) %>% 
  dplyr::summarize(snps = paste(snp_id, collapse = "€")) %>% 
  tidyr::separate(snps, c("snp1","snp2"), sep ="€")

#Calculate LD for each pair
snps = union(snps_per_gene$snp1, snps_per_gene$snp2)
selected_snps = vcf_file$genotypes[snps, ]
ld_mat = dplyr::group_by(snps_per_gene, gene_id) %>% 
  dplyr::mutate(ld = snpgdsLDpair(selected_snps[snp1,],selected_snps[snp2,], method = "dprime")[1])
dplyr::filter(ld_mat, ld*ld < 0.2)
