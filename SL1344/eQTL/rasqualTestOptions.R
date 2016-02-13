library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")

lib_size = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb.txt")
lib_size_gc = importRasqualTableOld("results/SL1344/rasqual/output/chr11_naive_100kb_gc.txt")
lib_size_gc_svd = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_svd.txt")
rle = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_RLE.txt")
rle_gc = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_RLE.txt")
lib_size_gc_PEER = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER.txt")
lib_size_gc_PEER_renamed = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_renamed.txt")
peer3 = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_n3.txt")

pop_only = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_n3_pop.txt")
pop_only_cov = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_pop.txt")
pop_only_svd = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_svd_pop.txt")


lib_size_min = findMinimalSnpPvalues(lib_size)
lib_size_gc_min = findMinimalSnpPvalues(lib_size_gc)
lib_size_gc_svd_min = findMinimalSnpPvalues(lib_size_gc_svd)
rle_min = findMinimalSnpPvalues(rle)
rle_gc_min = findMinimalSnpPvalues(rle_gc)
peer_gc_min = findMinimalSnpPvalues(lib_size_gc_PEER)
peer_gc_min_renamed = findMinimalSnpPvalues(lib_size_gc_PEER_renamed)
peer3_min = findMinimalSnpPvalues(peer3) 

pop_only_min = findMinimalSnpPvalues(pop_only)
pop_only_cov_min = findMinimalSnpPvalues(pop_only_cov)
pop_only_svd_min = findMinimalSnpPvalues(pop_only_svd)


#Explore results from eigenMT
eigen_res = eigenMTImportResults("results/SL1344/eigenMT/exampleOut.txt")


#Compare these results to Fastqtl
no_cov = importFastQTLTable("results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_perm.txt.gz") %>% findMinimalSnpPvalues()
cov = importFastQTLTable("results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_cov_perm.txt.gz") %>% findMinimalSnpPvalues()

#Peer covariate effect
a = dplyr::left_join(peer_gc_min_renamed, lib_size_gc_min, by = "gene_id")


#Peer vs fastqtl
b = dplyr::left_join(peer_gc_min_renamed, cov, by = "gene_id")
plot(-log(b$p_nominal.x,10), -log(b$p_nominal.y,10), xlim = c(0,70), ylim = c(0,70))


#External genotype matrix
eigen_res = eigenMTImportResults("results/SL1344/eigenMT/output/exampleOut.txt")
eigen_res_external = eigenMTImportResults("results/SL1344/eigenMT/output/test_out.txt")


#AcLDL dataset
lib_size_gc = importRasqualTable("results/acLDL/rasqual/output/chr11_naive_100kb_gc.txt") %>% findMinimalSnpPvalues() %>%
  dplyr::filter(p_fdr < 0.1)
lib_size_gc_PEER = importRasqualTable("results/acLDL/rasqual/output/chr11_naive_100kb_gc_PEER.txt") %>% findMinimalSnpPvalues() %>%
  dplyr::filter(p_fdr < 0.1)
lib_size_gc_svd = importRasqualTable("results/acLDL/rasqual/output/chr11_naive_100kb_gc_svd.txt") %>% findMinimalSnpPvalues() %>%
  dplyr::filter(p_fdr < 0.1)

