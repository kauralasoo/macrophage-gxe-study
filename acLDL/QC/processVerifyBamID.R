library("devtools")
library("dplyr")
load_all("macrophage-gxe-study/seqUtils/")

#Load verifyBamId results from disk
sample_names = read.table("macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_all.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
verify_bam_id = loadVerifyBamID(sample_names, "STAR/acLDL")

write.table(verify_bam_id, "macrophage-gxe-study/data/sample_lists/acLDL/acLDL_true_genotypes.txt")