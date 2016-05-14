library("GenomicRanges")

#Import credible sets
credible_sets = readRDS("results/SL1344/eQTLs/rasqual_credible_sets.rds")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Create a df of SNP positions
snp_pos_df = vcf_file$snpspos %>% 
  dplyr::transmute(seqnames = chr, start = pos, end = pos, strand = "+", snp_id = snpid)

#Convert credible sets to GRangesLists
#(This can take a long time)
cs_granges_list = purrr::map(credible_sets, 
                             ~purrr::map(.,~dplyr::filter(snp_pos_df, snp_id %in% .$snp_id) %>% 
                                           dataFrameToGRanges()) %>% 
                               GRangesList())
