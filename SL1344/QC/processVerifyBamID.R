library("dplyr")
library("devtools")
load_all("../seqUtils/")

#Load verifyBamID results from disk
sample_names = read.table("macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
verify_bam_id = loadVerifyBamID(sample_names, "STAR/SL1344")
write.table(verify_bam_id, "macrophage-gxe-study/data/sample_lists/SL1344/SL1344_true_genotype.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


