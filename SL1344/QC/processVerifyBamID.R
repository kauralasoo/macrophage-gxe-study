library("devtools")
load_all("macrophage-gxe-study/seqUtils/")

#Load verifyBamID results from disk
sample_names = read.table("fastq/SL1344_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
verify_bam_id = loadVerifyBamID(sample_names, "STAR/SL1344")
write.table(verify_bam_id, "macrophage-gxe-study/data/covariates/verifyBamID_true_genotype.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


