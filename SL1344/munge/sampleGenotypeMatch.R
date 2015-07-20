library("dplyr")
library("ggplot2")
library("tidyr")
library("devtools")
load_all("macrophage-gxe-study/seqUtils/")

line_meta = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")

#Import sample names
sample_names = read.table("sample_list.txt", stringsAsFactors = FALSE, comment = "")[,1]

#Match sample ids to genotypes
design = constructDesignMatrix_SL1344(sample_names)

gt_df = dplyr::select(line_meta, donor, line_id, genotype_id)
gt_map = dplyr::left_join(design, gt_df, by = "donor") %>% 
  dplyr::select(sample_id, genotype_id) %>% 
  unique()

write.table(gt_map, "genotypes/SL1344/SL1344_batch2_map.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)