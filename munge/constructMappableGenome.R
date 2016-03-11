
#Import blacklist and change chromosome names
blacklist = read.table("../../annotations/blacklists/wgEncodeDacMapabilityConsensusExcludable.hg38.bed")
chr_map = read.table("../macrophage-gxe-study/macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt")
colnames(chr_map) = c("V1", "V7")
blacklist_new_names = dplyr::filter(blacklist, V1 %in% chr_map$V1) %>%
  dplyr::left_join(chr_map, by = "V1") %>% 
  dplyr::select(-V1) %>%
  dplyr::select(V7, everything())
blacklist_gr = GRanges(seqnames = blacklist_new_names$V7, ranges = IRanges(start = blacklist_new_names$V2, end = blacklist_new_names$V3))

#Import gaps
gaps = read.table("../macrophage-gxe-study/macrophage-gxe-study/data/liftOver_genotypes/GRCh38_gaps.bed") %>%
  dplyr::filter(V1 %in% chr_map$V7)
gaps_gr = GRanges(seqnames = gaps$V1, ranges = IRanges(start = gaps$V2, end = gaps$V3))

#Import chromosome lengths
chr_lengths = read.table("../macrophage-gxe-study/macrophage-gxe-study/data/liftOver_genotypes/GRCh38_chromosome_lengths.txt") %>%
  dplyr::filter(V1 %in% chr_map$V7)
mappable_regions = GRanges(seqnames = chr_lengths$V1, ranges = IRanges(start = 1, end = chr_lengths$V2))

#Construct list of mappable regions
filtered_genome = setdiff(setdiff(mappable_regions, gaps_gr), blacklist_gr)


