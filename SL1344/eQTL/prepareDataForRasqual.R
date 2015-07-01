library("dplyr")
library("devtools")
library("GenomicFeatures")
load_all("macrophage-gxe-study/seqUtils/")

#Import data
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
txdb79 = loadDb("../../annotations/GRCh38/genes/TranscriptDb_GRCh38_79.db")
gencode_basic_transcript = read.table("../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.txt", 
                                      header = TRUE, stringsAsFactors = FALSE) %>% tbl_df()

#Make the design matrix
design = dplyr::filter(expression_dataset$design, !(donor == "ffdp")) %>% #ffdp identical to fpdl
  dplyr::filter(!(donor == "fpdj" & replicate == 2)) #fpdj has two replicates
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#Make a table of genotype ids
genotype_ids = dplyr::filter(sample_meta, condition == "A") %>% dplyr::select(line_id, genotype_id)
write.table(genotype_ids,"genotypes/SL1344/SL1344_genotype_list.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

#Make a table of sample-genotype match
sg_match = dplyr::filter(sample_meta) %>% dplyr::select(sample_id, genotype_id)
write.table(sg_match,"genotypes/SL1344/SL1344_sample_genotype_map.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

#Extract exon start and end coordinates from txdb 
exons = exonsBy(txdb79, by = "tx", use.names = TRUE)
gene_ids = unique(sort(gencode_basic_transcript$ensembl_gene_id))
gene_id_list = as.list(gene_ids)
names(gene_id_list) = gene_ids

gene_exons = lapply(gene_id_list, findGeneExons, gencode_basic_transcript, exons)
