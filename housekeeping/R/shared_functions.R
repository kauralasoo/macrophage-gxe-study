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

#Return a list of lists of paths to QTL mapping results
qtlResults <- function(){
  atac_rasqual = list(naive = "/Users/alasoo/databases/ATAC/rasqual/naive_100kb.sorted.txt.gz",
                      IFNg = "/Users/alasoo/databases/ATAC/rasqual/IFNg_100kb.sorted.txt.gz",
                      SL1344 = "/Users/alasoo/databases/ATAC/rasqual/SL1344_100kb.sorted.txt.gz",
                      IFNg_SL1344 = "/Users/alasoo/databases/ATAC/rasqual/IFNg_SL1344_100kb.sorted.txt.gz")
  
  atac_fastqtl = list(naive = "/Users/alasoo/databases/ATAC/fastqtl/naive_500kb_pvalues.sorted.txt.gz",
                      IFNg = "/Users/alasoo/databases/ATAC/fastqtl/IFNg_500kb_pvalues.sorted.txt.gz",
                      SL1344 = "/Users/alasoo/databases/ATAC/fastqtl/SL1344_500kb_pvalues.sorted.txt.gz",
                      IFNg_SL1344 = "/Users/alasoo/databases/ATAC/fastqtl/IFNg_SL1344_500kb_pvalues.sorted.txt.gz")
  
  rna_rasqual = list(naive = "/Users/alasoo/databases/SL1344/rasqual/naive_500kb.sorted.txt.gz",
                     IFNg = "/Users/alasoo/databases/SL1344/rasqual/IFNg_500kb.sorted.txt.gz",
                     SL1344 = "/Users/alasoo/databases/SL1344/rasqual/SL1344_500kb.sorted.txt.gz",
                     IFNg_SL1344 = "/Users/alasoo/databases/SL1344/rasqual/IFNg_SL1344_500kb.sorted.txt.gz")
  
  rna_fastqtl = list(naive = "/Users/alasoo/databases/SL1344/fastqtl/naive_500kb_pvalues.sorted.txt.gz",
                     IFNg = "/Users/alasoo/databases/SL1344/fastqtl/IFNg_500kb_pvalues.sorted.txt.gz",
                     SL1344 = "/Users/alasoo/databases/SL1344/fastqtl/SL1344_500kb_pvalues.sorted.txt.gz",
                     IFNg_SL1344 = "/Users/alasoo/databases/SL1344/fastqtl/IFNg_SL1344_500kb_pvalues.sorted.txt.gz") 
  
  result = list(atac_fastqtl = atac_fastqtl, atac_rasqual = atac_rasqual,
                rna_fastqtl = rna_fastqtl, rna_rasqual = rna_rasqual)
  return(result)
}

figureNames<- function(){
  data_frame(condition_name = c("naive", "IFNg", "SL1344", "IFNg_SL1344"), 
             figure_name = factor(c("N","I","S","I+S"), levels = c("N","I","S","I+S")))
  }
