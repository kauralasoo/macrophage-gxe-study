library("plyr")
library("dplyr")
library("MatrixEQTL")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

#Extract all SNP-gene pairs
min_pvalues_list = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list[1:2], function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalues_df = ldply(min_pvalues_hits, .id = "condition_name")
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()

#Load data
vcf_file = readRDS("../macrophage-gxe-study/results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata = dplyr::mutate(atac_list$sample_metadata, genotype_id = ifelse(genotype_id == "HPSI1213i-nusw_2", "HPSI1213i-nusw_1", genotype_id)) #Fix this temportary bug in genotypes vcf file

#Filter SNPs based on R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
joint_pairs_filtered = filterHitsR2(joint_pairs, genotypes, 0.8)

#Test pairwise interactions using matrixeQTL
interaction_df = matrixeqtlTestInteraction(c("naive", "IFNg"), joint_pairs_filtered, atac_list, vcf_file)

#Find effect sizes for these genes from RASQUAL data
selected_pvalue_list = readRDS("results/ATAC/QTLs/rasqual_selected_pvalues.rds")
interaction_effects = filterInteractionResults(c("naive", "IFNg"), interaction_df, selected_pvalue_list, fdr_thresh = 0.1)
interaction_effects_filtered = dplyr::filter(interaction_effects, abs_beta_min < 1, abs(beta_diff) > 1)




#Make database connections
naive_sql <- src_sqlite("~/projects/macrophage-chromatin/naive_50kb.sqlite3") %>% tbl("rasqual")
ifng_sql <- src_sqlite("~/projects/macrophage-chromatin/IFNg_50kb.sqlite3") %>% tbl("rasqual")

naive_tbl = dplyr::select(tbl_df(naive_res), gene_id, snp_id, beta, p_nominal) %>% dplyr::rename(beta_A = beta, p_nominal_A = p_nominal)
ifng_tbl = dplyr::select(tbl_df(ifng_res), gene_id, snp_id, beta, p_nominal) %>% dplyr::rename(beta_B = beta, p_nominal_B = p_nominal)

res = dplyr::left_join(interaction_candidates, naive_tbl, by = c("gene_id", "snp_id")) %>% 
  tbl_df() %>% 
  dplyr::left_join(ifng_tbl, by =c("gene_id","snp_id")) %>% 
  dplyr::mutate(beta_diff = beta_B - beta_A) %>%
  dplyr::mutate(abs_beta_min = pmin(abs(beta_A), abs(beta_B))) %>%
  dplyr::group_by(gene_id, snp_id) %>%
  dplyr::mutate(min_condition = which.min(c(abs(beta_A), abs(beta_B))))

filter

#Test overlap between credible sets of associated SNPs
gene_id_list = idVectorToList(different_lead_snps$gene_id)
overlaps = lapply(gene_id_list, testCredibleSetOverlapConditions, naive_sql, ifng_sql,27,27)
cs_overlap_df = plyr::ldply(overlaps, .id = "gene_id") %>%
  dplyr::rename(overlap = V1)

#Test for colocalisation between different lead SNPs
gene_id_list = idVectorToList(joint_hits$gene_id)
coloc_res = lapply(gene_id_list, testColocConditons, naive_sql, ifng_sql, 27, 27)
coloc_df = plyr::ldply(coloc_res, .id = "gene_id")
coloc_df$model = apply(coloc_df[,3:7], 1,which.max)



#Examples
naive_pvals = fetchSQLite(naive_sql, selected_gene_id = "ATAC_peak_280095", selected_snp_id = "rs17198394")
ifng_pvals = fetchSQLite(ifng_sql, selected_gene_id = "ATAC_peak_280095", selected_snp_id = "rs17198394")

col = testColoc(naive_pvals, ifng_pvals, 27, 27)

naive_post = addAssociationPosterior(naive_pvals, 27)
ifng_post = addAssociationPosterior(naive_pvals, 27)

naive_pvals = fetchSQLite(naive_sql, selected_gene_id = "ATAC_peak_280095") %>% addAssociationPosterior(27)
ifng_pvals = fetchSQLite(ifng_sql, selected_gene_id = "ATAC_peak_280095") %>% addAssociationPosterior(27)
plot(naive_pvals$pos, naive_pvals$posterior)
plot(ifng_pvals$pos, ifng_pvals$posterior)


#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
#Import genotypes
vcf_file = readRDS("../macrophage-gxe-study/results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")
atac_metadata = readRDS("macrophage-chromatin/data/SL1344/compiled_atac_metadata.rds") %>% 
  dplyr::semi_join(atac_list$design, by = "sample_id")

ctsc_enhancer = plotEQTL("ATAC_peak_158153", "rs12876177", log(atac_list$tpm + 0.1,2), vcf_file$genotypes, 
                         atac_list$sample_metadata, atac_list$gene_metadata)

plotEQTL("ATAC_peak_252666", "rs10872076", atac_list$cqn, vcf_file$genotypes, 
         atac_list$sample_metadata, atac_list$gene_metadata)
