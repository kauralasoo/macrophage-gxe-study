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
            "genotypes/SL1344/SL1344_sg_map_A.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(dplyr::filter(sg_match, condition == "B") %>% dplyr::select(-condition),
            "genotypes/SL1344/SL1344_sg_map_B.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(dplyr::filter(sg_match, condition == "C") %>% dplyr::select(-condition),
            "genotypes/SL1344/SL1344_sg_map_C.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(dplyr::filter(sg_match, condition == "D") %>% dplyr::select(-condition),
            "genotypes/SL1344/SL1344_sg_map_D.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
