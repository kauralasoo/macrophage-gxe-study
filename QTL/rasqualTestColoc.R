library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")

#Import minimal p-values for each peak
min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")

#Load minimal SNPs from disk
naive_hits = dplyr::filter(min_pvalues$naive, p_fdr < 0.1) %>% dplyr::transmute(gene_id, naive_snp_id = snp_id)
ifng_hits = dplyr::filter(min_pvalues$IFNg, p_fdr < 0.1) %>% dplyr::transmute(gene_id, IFNg_snp_id = snp_id)

joint_hits = dplyr::full_join(naive_hits, ifng_hits, by = "gene_id")
different_lead_snps = dplyr::filter(joint_hits, naive_snp_id != IFNg_snp_id)

#Test for colocalization for different lead SNPs
naive_sql <- src_sqlite("~/projects/macrophage-chromatin/naive_50kb.sqlite3") %>% tbl("rasqual")
ifng_sql <- src_sqlite("~/projects/macrophage-chromatin/IFNg_50kb.sqlite3") %>% tbl("rasqual")

naive_pvals = fetchSQLite(naive_sql, selected_gene_id = "ATAC_peak_55686")
ifng_pvals = fetchSQLite(ifng_sql, selected_gene_id = "ATAC_peak_55686")

col = testColoc(naive_pvals, ifng_pvals, 27, 27)
naive_post = addAssociationPosterior(naive_pvals, 27)

testColocConditons <- function(gene_id, db_table1, db_table2, n1, n2, p1, p2, p12){
  p_values_1 = fetchSQLite(db_table1, selected_gene_id = gene_id)
  p_values_2 = fetchSQLite(db_table2, selected_gene_id = gene_id)
  
  coloc = testColoc(p_values_1, p_values_2, n1, n2, p1, p2, p12)
  return(coloc$summary)
}

a = testColocConditons("ATAC_peak_3705", naive_sql, ifng_sql, 27, 27)

#Test for colocalisation between different lead SNPs
gene_id_list = idVectorToList(joint_hits$gene_id)
coloc_res = lapply(gene_id_list, testColocConditons, naive_sql, ifng_sql, 27, 27, 1e-3,1e-3, 5e-4)
coloc_df = plyr::ldply(coloc_res, .id = "gene_id")
coloc_df$model = apply(coloc_df[,3:7], 1,which.max)

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
#Import genotypes
vcf_file = readRDS("../macrophage-gxe-study/results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")
atac_metadata = readRDS("macrophage-chromatin/data/SL1344/compiled_atac_metadata.rds") %>% 
  dplyr::semi_join(atac_list$design, by = "sample_id")

ctsc_enhancer = plotEQTL("ATAC_peak_245806", "rs4714033", atac_list$exprs_cqn, vcf_file$genotypes, 
                         atac_metadata, atac_list$gene_metadata)
