library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")

lib_size = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb.txt")
lib_size_gc = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc.txt")
lib_size_gc_svd = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_svd.txt")
rle = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_RLE.txt")
rle_gc = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_RLE.txt")
lib_size_gc_PEER = importRasqualTable("results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER.txt")

lib_size_min = findMinimalSnpPvalues(lib_size)
lib_size_gc_min = findMinimalSnpPvalues(lib_size_gc)
lib_size_gc_svd_min = findMinimalSnpPvalues(lib_size_gc_svd)
rle_min = findMinimalSnpPvalues(rle)
rle_gc_min = findMinimalSnpPvalues(rle_gc)
peer_gc_min = findMinimalSnpPvalues(lib_size_gc_PEER)

#Explore results from eigenMT
eigen_res = eigenMTImportResults("results/SL1344/eigenMT/exampleOut.txt")
