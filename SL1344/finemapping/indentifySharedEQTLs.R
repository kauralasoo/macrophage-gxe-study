library("GenomicRanges")
load_all("../seqUtils/")

#Import credible sets
credible_sets = readRDS("results/SL1344/eQTLs/rasqual_credible_sets.rds")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Create a df of SNP positions
snp_pos_df = vcf_file$snpspos %>% 
  dplyr::transmute(seqnames = chr, start = pos, end = pos, strand = "+", snp_id = snpid)

#Convert credible sets into gigantic data frame
credible_sets_df = purrr::map_df(credible_sets, ~purrr::map_df(., 
                                 ~dplyr::mutate(.,chr = as.character(chr))) %>% 
                                 dplyr::filter(chr != "X"), .id = "condition_name")

#Construct granges object
credible_sets_granges = credible_sets_df %>% 
  dplyr::transmute(condition_name, gene_id, snp_id, chr, seqnames = chr, start = pos, end = pos, strand = "+") %>% 
  dataFrameToGRanges()

#Find overlaps between credible sets
olaps = findOverlaps(credible_sets_granges)
queries = credible_sets_granges[queryHits(olaps),] %>% 
  elementMetadata() %>% as.data.frame() %>% 
  dplyr::transmute(master_condition = condition_name, master_id = gene_id, chr1 = chr) %>% 
  tbl_df()
subjects = credible_sets_granges[subjectHits(olaps),] %>% 
  elementMetadata() %>% as.data.frame() %>% 
  dplyr::transmute(dependent_condition = condition_name, dependent_id = gene_id, chr2 = chr) %>% 
  tbl_df()

#Identify all genes that share at least one SNP in their credible set
pairwise_shared = dplyr::bind_cols(queries, subjects) %>% 
  unique() %>% 
  dplyr::filter(chr1 == chr2) %>% 
  dplyr::filter(master_id != dependent_id) %>% 
  dplyr::select(master_id, master_condition, dependent_id, dependent_condition) %>% 
  unique()

a = credible_sets$naive$ENSG00000164308$snp_id
b = credible_sets$IFNg$ENSG00000113441$snp_id
length(intersect(a,b))/length(union(a,b))

#Convert shared peaks into clusters
clusters = constructClustersFromGenePairs(pairwise_shared)

#Use Coloc to test if the QTLs are indeed shared

#Test for GxGxE interactions





