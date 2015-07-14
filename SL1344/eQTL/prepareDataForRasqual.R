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
#Correct sample swap between vorx and zuta
line_metadata$line_id[c(which(line_metadata$line_id == "zuta_1"), which(line_metadata$line_id == "vorx_1"))] = c("vorx_1", "zuta_1")
line_metadata$donor[c(which(line_metadata$donor == "zuta"), which(line_metadata$donor == "vorx"))] = c("vorx", "zuta")

#Make the design matrix
design = dplyr::filter(expression_dataset$design, !(donor == "ffdp")) %>% #ffdp identical to fpdl
  dplyr::filter(!(donor == "fpdj" & replicate == 2)) #fpdj has two replicates
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#Make a table of genotype ids
genotype_ids = dplyr::filter(sample_meta, condition == "A") %>% dplyr::select(line_id, genotype_id)
write.table(genotype_ids,"genotypes/SL1344/SL1344_genotype_list.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

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

#Export read counts in each condtion
cond_A_counts = expression_dataset$exprs_counts[,dplyr::filter(sg_match, condition == "A")$sample_id]
write.table(cond_A_counts, "rasqual/input/cond_A_counts.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
cond_B_counts = expression_dataset$exprs_counts[,dplyr::filter(sg_match, condition == "B")$sample_id]
write.table(cond_B_counts, "rasqual/input/cond_B_counts.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
cond_C_counts = expression_dataset$exprs_counts[,dplyr::filter(sg_match, condition == "C")$sample_id]
write.table(cond_C_counts, "rasqual/input/cond_C_counts.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
cond_D_counts = expression_dataset$exprs_counts[,dplyr::filter(sg_match, condition == "D")$sample_id]
write.table(cond_D_counts, "rasqual/input/cond_D_counts.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)

#Estimate size factors for each sample
cond_A_factors = calculateNormFactors(cond_A_counts, method = "RLE")
write.table(cond_A_factors, "rasqual/input/cond_A_factors.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
cond_B_factors = calculateNormFactors(cond_B_counts, method = "RLE")
write.table(cond_B_factors, "rasqual/input/cond_B_factors.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
cond_C_factors = calculateNormFactors(cond_C_counts, method = "RLE")
write.table(cond_C_factors, "rasqual/input/cond_C_factors.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)
cond_D_factors = calculateNormFactors(cond_D_counts, method = "RLE")
write.table(cond_D_factors, "rasqual/input/cond_D_factors.txt", row.names = TRUE, col.names = FALSE, sep ="\t", quote = FALSE)


