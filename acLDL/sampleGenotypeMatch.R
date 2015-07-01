library("dplyr")
library("ggplot2")
library("tidyr")
library("devtools")
load_all("macrophage-gxe-study/seqUtils/")

#Import Compiled metadata
line_meta = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")
#Import acLDL sample names
acldl_names = read.table("fastq/acLDL_names_all.txt", stringsAsFactors = FALSE, comment = "")[,1]

#Match sample ids to genotypes
acldl_design = constructDesignMatrix_acLDL(acldl_names)
gt_df = dplyr::select(line_meta, donor, line_id, genotype_id)
gt_map = dplyr::left_join(acldl_design, gt_df, by = "donor") %>% 
  dplyr::select(sample_id, genotype_id) %>% 
  unique()

#Save map to disk
write.table(gt_map, "genotypes/acLDL/acLDL_sample_genotype_map.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)