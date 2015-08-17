library("dplyr")
library("devtools")
load_all("../macrophage-gxe-study/macrophage-gxe-study/seqUtils/")

#Load verifyBamID results from disk
sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
verify_bam_id = loadVerifyBamID(sample_names, "processed/SL1344/")
write.table(verify_bam_id, "results/ATAC/QC/ATAC_sample_genotype_match.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


