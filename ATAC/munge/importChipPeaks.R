library("devtools")
library("dplyr")
library("rtracklayer")
load_all("../seqUtils/")

#### Import peaks for enrichment in naive condition ####
#Load Rehli dataset
sample_names = c("MAC_PU1","MAC_CEBPbeta","CTCF_MAC")
rehli_peaks = loadNarrowPeaks("results/Rehli/peak_calls/", sample_names, sub_dir = FALSE)
elementMetadata(rehli_peaks$MAC_PU1)$name = "PU.1_Pham"
elementMetadata(rehli_peaks$MAC_CEBPbeta)$name = "CEBPb_Pham"
elementMetadata(rehli_peaks$CTCF_MAC)$name = "CTCF_Pham"
rehli_peaks = purrr::map(rehli_peaks, function(x){elementMetadata(x) = elementMetadata(x)[,"name",drop=FALSE]; return(x)})

#Load OCallaghan datset
cebpb_ctrl_names = c("CEBPbeta_ctrl_201","CEBPbeta_ctrl_202","CEBPbeta_ctrl_203","CEBPbeta_ctrl_204")
CEBPbeta_ctrl_peaks = loadNarrowPeaks("results/OCallaghan/peak_calls/", cebpb_ctrl_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
elementMetadata(CEBPbeta_ctrl_peaks)$name = "CEBPb_Reschen"

#Load Schultze dataset
schultze_names = c("PU1_naive")
schultze_peaks = loadNarrowPeaks("results/Schultze/peak_calls/", schultze_names, sub_dir = FALSE)
PU1_Schultze = schultze_peaks$PU1_naive
PU1_Schultze = PU1_Schultze[PU1_Schultze$score > 100,]
elementMetadata(PU1_Schultze)$name = "PU.1_Schmidt"
elementMetadata(PU1_Schultze) = elementMetadata(PU1_Schultze)[,"name", drop=FALSE]

#Join peaks together
joint_peaks = c(rehli_peaks$MAC_PU1, rehli_peaks$MAC_CEBPbeta, rehli_peaks$CTCF_MAC, CEBPbeta_ctrl_peaks, PU1_Schultze)
rtracklayer::export.bed(joint_peaks, "results/ATAC/ChIP_enrichment/naive_combined_peaks.bed")
saveRDS(joint_peaks, "results/ATAC/ChIP_enrichment/naive_combined_peaks.rds")


#### Load Ivashkiv peaks ####
sample_names = c("STAT1_rep1_A","STAT1_rep1_B","STAT1_rep1_C","STAT1_rep1_D","STAT1_rep2_B","STAT1_rep2_D",
                 "IRF1_A","IRF1_B","IRF1_E","IRF1_F")
ivashkiv_peaks = loadNarrowPeaks("results/Ivashkiv/peak_calls/", sample_names, sub_dir = FALSE)
saveRDS(ivashkiv_peaks, "results/ATAC/ChIP_enrichment/ivashkiv_peaks.rds")


#### Load and process peaks from the Wong et al dataset ####

#naive CIITA
ciita_naive_names = c("MO_rep1_naive_CIITA","MO_rep2_naive_CIITA")
ciita_naive_peaks = loadNarrowPeaks("results/Knight/peak_calls/", ciita_naive_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
elementMetadata(ciita_naive_peaks)$name = "CIITA_naive"

#IFNg CIITA
ciita_ifng_names = c("MO_rep1_IFNg_CIITA","MO_rep2_IFNg_CIITA")
ciita_ifng_peaks = loadNarrowPeaks("results/Knight/peak_calls/", ciita_ifng_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
elementMetadata(ciita_ifng_peaks)$name = "CIITA_IFNg"

#naive RFX5
rfx5_naive_names = c("MO_rep1_naive_RFX5","MO_rep2_naive_RFX5")
rfx5_naive_peaks = loadNarrowPeaks("results/Knight/peak_calls/", rfx5_naive_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
elementMetadata(rfx5_naive_peaks)$name = "RFX5_naive"


rfx5_ifng_names = c("MO_rep1_IFNg_RFX5","MO_rep2_IFNg_RFX5")
rfx5_ifng_peaks = loadNarrowPeaks("results/Knight/peak_calls/", rfx5_ifng_names, sub_dir = FALSE) %>% 
  filterOverlaps(minOverlapCount = 2) %>% listUnion()
elementMetadata(rfx5_ifng_peaks)$name = "RFX5_IFNg"

joint_peaks = c(ciita_naive_peaks, ciita_ifng_peaks, rfx5_naive_peaks, rfx5_ifng_peaks)
export.bed(joint_peaks, "results/ATAC/ChIP_enrichment/CIITA-RFX5_joint_peaks.bed")
saveRDS(joint_peaks, "results/ATAC/ChIP_enrichment/CIITA-RFX5_joint_peaks.rds")


#### Make consensus peaks ####
chr_name = read.table("../macrophage-gxe-study/macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt", stringsAsFactors = FALSE)

#IRF1
irf_names = c("IRF1_A","IRF1_B","IRF1_E","IRF1_F")
irf_peaks = loadNarrowPeaks("processed/Ivashkiv/", irf_names)
irf_peaks = makeUnionPeaks(irf_peaks, chr_name[,1], "IRF1_peak_")
rtracklayer::export.gff3(irf_peaks, "annotations/IRF1_joint_peaks.gff3")

#STAT1
stat1_names = c("STAT1_rep1_A","STAT1_rep1_B","STAT1_rep1_C","STAT1_rep1_D")
stat1_peaks = loadNarrowPeaks("processed/Ivashkiv/", stat1_names)
stat1_peaks = makeUnionPeaks(stat1_peaks, chr_name[,1], "STAT1_peak_")
rtracklayer::export.gff3(stat1_peaks, "annotations/STAT1_joint_peaks.gff3")
