library("plyr")
library("dplyr")
library("tidyr")
library("devtools")
library("GenomicFeatures")
load_all("macrophage-gxe-study/seqUtils/")

#Import data
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") %>% #Line metadata
  dplyr::filter(status == "Success")

#Make the design matrix
design = dplyr::filter(expression_dataset$design, !(donor == "fpdj")) %>% tbl_df() %>% #Remove all fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) #Remove second ougl sample (dium)
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#Match samples to genotypes in each condition separately
sg_match = dplyr::filter(sample_meta) %>% dplyr::select(genotype_id, sample_id, condition)

write.table(dplyr::filter(sg_match, condition == "A") %>% dplyr::select(-condition),
            "rasqual/input/SL1344_sg_map_A.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(dplyr::filter(sg_match, condition == "B") %>% dplyr::select(-condition),
            "rasqual/input/SL1344_sg_map_B.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(dplyr::filter(sg_match, condition == "C") %>% dplyr::select(-condition),
            "rasqual/input/SL1344_sg_map_C.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(dplyr::filter(sg_match, condition == "D") %>% dplyr::select(-condition),
            "rasqual/input/SL1344_sg_map_D.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

#Export read counts in each condtion (Including binary format)
cond_A_counts = expression_dataset$exprs_counts[,dplyr::filter(sg_match, condition == "A")$sample_id]
write.table(cond_A_counts, "rasqual/input/counts/cond_A_counts.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
writeBin(as.double(c(t(cond_A_counts))), "rasqual/input/counts/cond_A_counts.bin")

cond_B_counts = expression_dataset$exprs_counts[,dplyr::filter(sg_match, condition == "B")$sample_id]
write.table(cond_B_counts, "rasqual/input/counts/cond_B_counts.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
writeBin(as.double(c(t(cond_B_counts))), "rasqual/input/counts/cond_B_counts.bin")

cond_C_counts = expression_dataset$exprs_counts[,dplyr::filter(sg_match, condition == "C")$sample_id]
write.table(cond_C_counts, "rasqual/input/counts/cond_C_counts.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
writeBin(as.double(c(t(cond_C_counts))), "rasqual/input/counts/cond_C_counts.bin")

cond_D_counts = expression_dataset$exprs_counts[,dplyr::filter(sg_match, condition == "D")$sample_id]
write.table(cond_D_counts, "rasqual/input/counts/cond_D_counts.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
writeBin(as.double(c(t(cond_D_counts))), "rasqual/input/counts/cond_D_counts.bin")

#Estimate size factors for each sample
##TODO: Implement GC-content correction into size factor caluclation
cond_A_factors = calculateNormFactors(cond_A_counts, method = "RLE")
write.table(cond_A_factors, "rasqual/input/counts/cond_A_factors.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
writeBin(as.double(c(t(cond_A_factors))), "rasqual/input/counts/cond_A_factors.bin")

cond_B_factors = calculateNormFactors(cond_B_counts, method = "RLE")
write.table(cond_B_factors, "rasqual/input/counts/cond_B_factors.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
writeBin(as.double(c(t(cond_B_factors))), "rasqual/input/counts/cond_B_factors.bin")

cond_C_factors = calculateNormFactors(cond_C_counts, method = "RLE")
write.table(cond_C_factors, "rasqual/input/counts/cond_C_factors.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
writeBin(as.double(c(t(cond_C_factors))), "rasqual/input/counts/cond_C_factors.bin")

cond_D_factors = calculateNormFactors(cond_D_counts, method = "RLE")
write.table(cond_D_factors, "rasqual/input/counts/cond_D_factors.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
writeBin(as.double(c(t(cond_D_factors))), "rasqual/input/counts/cond_D_factors.bin")


