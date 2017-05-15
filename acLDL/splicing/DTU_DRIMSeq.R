library("dplyr")
library("DRIMSeq")
library("devtools")
library("optparse")
library("SummarizedExperiment")
load_all("../seqUtils/")

##Parse command-line options
option_list <- list(
  make_option(c("-p", "--phenotype"), type="character", default=NULL,
              help="Type of transcript quantification used for DTU testing.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Import SummarizedExperiments for all phenotypes
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
se_revised = readRDS("results/acLDL/acLDL_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/acLDL/acLDL_leafcutter_counts.rds")
se_list = list(ensembl_87 = se_ensembl, revisedAnnotation = se_revised, leafcutter = se_leafcutter)

#Select one phenotype
phenotype = opt$p
selected_se = se_list[[phenotype]]
event_metadata = rowData(selected_se) %>% tbl_df2()
unique_genes = unique(event_metadata$gene_id)

####### Get batch string from stdin ######
f <- file("stdin")
open(f)
batch_string = readLines(f)
close(f)
####### END #######

#Pefromed DTU test on LeafCutter counts
#Construct splicing dataset
counts_matrix = assays(selected_se)$counts %>% tbl_df2() %>%
  dplyr::mutate(gene_id = rowData(selected_se)$gene_id, feature_id = rowData(selected_se)$transcript_id) %>%
  dplyr::select(gene_id, feature_id, everything())
sample_matrix = colData(selected_se) %>% 
  tbl_df2() %>% 
  dplyr::transmute(sample_id, group = condition_name)

#Create DRIMSeq object
splicing_data <- DRIMSeq::dmDSdata(counts = as.data.frame(counts_matrix), 
                                   samples = as.data.frame(sample_matrix))

#Filter
splicing_data <- DRIMSeq::dmFilter(splicing_data, min_samps_gene_expr = 70, min_samps_feature_expr = 70,
                          min_samps_feature_prop = 0)

#Select genes from a single batch
batch_genes = constructIdBatches(batch_string, unique_genes)
splicing_data_batch <- splicing_data[names(splicing_data) %in% batch_genes, ]

#Construct model matrix
design_full <- model.matrix(~group, data = DRIMSeq::samples(splicing_data_batch))

#Estimate precision and fit model
model_fit <- dmPrecision(splicing_data_batch, design = design_full, verbose = 1)
model_fit <- dmFit(model_fit, design = design_full, verbose = 1)

#Test DTU
design_null <- model.matrix(~ 1, data = DRIMSeq::samples(splicing_data_batch))
model_fit <- dmTest(model_fit, design = design_null)

#Save output from each batch
output_file = file.path("results/acLDL/diff_splicing", paste0("differential_events_batch",gsub(" ","_",batch_string), ".rds"))
saveRDS(model_fit, output_file)

