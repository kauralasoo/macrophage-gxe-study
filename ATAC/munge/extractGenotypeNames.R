library("devtools")
library("dplyr")
load_all("macrophage-chromatin/housekeeping/")

sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
sample_meta = readRDS("../macrophage-gxe-study/macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")

#Construct design matrix and add metadata
atac_design = constructDesignMatrix_ATAC(sample_names)
atac_meta = dplyr::left_join(atac_design, sample_meta, by = "donor")

#Save genotype names to disk
atac_genotypes = unique(atac_meta$genotype_id)
write.table(atac_genotypes, "macrophage-chromatin/data/SL1344/ATAC_gt_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
