library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")

#Import minimal p-values for each peak
min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")

#Load minimal SNPs from disk
naive_hits = dplyr::filter(min_pvalues$naive, p_fdr < 0.1) %>% dplyr::transmute(gene_id, naive_snp_id = snp_id)
ifng_hits = dplyr::filter(min_pvalues$IFNg, p_fdr < 0.1) %>% dplyr::transmute(gene_id, IFNg_snp_id = snp_id)

joint_pairs = rbind(dplyr::filter(min_pvalues$naive, p_fdr < 0.1), dplyr::filter(min_pvalues$IFNg, p_fdr < 0.1)) %>% 
  dplyr::select(gene_id, snp_id) %>% unique()

joint_hits = dplyr::full_join(naive_hits, ifng_hits, by = "gene_id")
different_lead_snps = dplyr::filter(joint_hits, naive_snp_id != IFNg_snp_id)

#Make database connections
naive_sql <- src_sqlite("~/projects/macrophage-chromatin/naive_50kb.sqlite3") %>% tbl("rasqual")
ifng_sql <- src_sqlite("~/projects/macrophage-chromatin/IFNg_50kb.sqlite3") %>% tbl("rasqual")

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

#Test pairwise interactions using matrixeQTL
vcf_file = readRDS("../macrophage-gxe-study/results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata = dplyr::mutate(atac_list$sample_metadata, genotype_id = ifelse(genotype_id == "HPSI1213i-nusw_2", "HPSI1213i-nusw_1", genotype_id)) #Fix this temportary bug in genotypes vcf file

design_matrix = dplyr::filter(atac_list$sample_metadata, condition_name %in% c("naive","IFNg"))
exp_matrix = atac_list$cqn[joint_hits$gene_id,design_matrix$sample_id]

snp_ids = union(naive_hits$naive_snp_id, ifng_hits$IFNg_snp_id)
geno_matrix = extractSubset(design_matrix, vcf_file$genotypes[snp_ids,], old_column_names = "genotype_id", new_column_names = "sample_id")
snpspos = dplyr::filter(vcf_file$snpspos, snpid %in% snp_ids) %>% as.data.frame()
genepos = constructMatrixEQTLGenePos(atac_list$gene_metadata) %>% dplyr::filter(geneid %in% joint_hits$gene_id)

#Construct condition covariate
cov = dplyr::mutate(design_matrix, condition_cov = ifelse(condition_name == "naive", 0, 1))$condition_cov
cov_matrix = t(as.matrix(cov))
colnames(cov_matrix) = design_matrix$sample_id

interaction_res = runMatrixEQTL(exp_matrix, geno_matrix, snpspos, genepos, cov_matrix, pvOutputThreshold = 1, model = modelLINEAR_CROSS)
interaction_table = interaction_res$cis$eqtls %>% dplyr::rename(gene_id = gene, snp_id = snps) %>% 
  dplyr::semi_join(joint_pairs, by = c("gene_id", "snp_id")) %>%
  dplyr::mutate(p_fdr = p.adjust(pvalue, "fdr")) %>%
  dplyr::arrange(p_fdr)
  

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

plotEQTL("ATAC_peak_158153", "rs13031964", atac_list$cqn, vcf_file$genotypes, 
         atac_list$sample_metadata, atac_list$gene_metadata)
