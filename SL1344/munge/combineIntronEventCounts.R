library("devtools")
library("dplyr")
library("readr")
load_all("macrophage-gxe-study/seqUtils/")

sample_names = read.table("macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata

condA_names = sample_names[grepl("_A", sample_names)]
design_matrix = constructDesignMatrix_SL1344(condA_names)
design = dplyr::filter(design_matrix, !(donor == "fpdj")) %>% tbl_df() %>% #Remove all fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) #Remove second ougl sample (dium)

sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#filter samples by condition
condA_counts = loadIntronEventCounts("STAR/SL1344/", design$sample_id, sub_dir = TRUE)
new_cols = c(colnames(condA_counts)[1:2], sample_meta$genotype_id)
colnames(condA_counts) = new_cols
write.table(condA_counts, "results/SL1344/condition_A_intron_event_counts.txt", sep ="\t", row.names = FALSE, quote = FALSE)