library("dplyr")

results1 = read.table(gzfile("eqtlbma/output/out_eqtlbma_avg_bfs.txt.gz"), sep = "\t", skip = 1, header = TRUE)
significant_genes = dplyr::filter(results1, gene.post > 0.9) %>% tbl_df() %>%
  dplyr::select(gene, snp, gene.post, snp.post.the, snp.post.an, best.config, post.best.config) %>%
  dplyr::group_by(gene) %>%
  dplyr::arrange(-snp.post.an) %>%
  dplyr::filter(row_number() == 1)
table(significant_genes$best.config)