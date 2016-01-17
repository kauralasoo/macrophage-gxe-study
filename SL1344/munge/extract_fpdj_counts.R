library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Load the raw eQTL dataset
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
fpdj_design = expression_dataset$design %>% dplyr::filter(donor %in% c("fpdj","nibo"))
fpdj_counts = expression_dataset$exprs_counts[,fpdj_design$sample_id]

write.table(fpdj_counts, "results/SL1344/fpdj_counts/count_matrix.txt", sep = "\t", quote = FALSE)
write.table(fpdj_design, "results/SL1344/fpdj_counts/design_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
