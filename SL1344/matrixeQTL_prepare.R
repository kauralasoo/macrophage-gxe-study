library("dplyr")
library("tidyr")
library("devtools")
library("cqn")
load_all("macrophage-gxe-study/seqUtils/")

#### Load raw counts from disk ####
data = read.table("results//SL1344//combined_counts.txt", stringsAsFactors = FALSE) #Read counts
#There seems to be a sample swap between coio_C and coio_D, fix that
indexes = c(which(colnames(data) == "coio_C"), which(colnames(data) == "coio_D"))
colnames(data)[indexes] = c("coio_D", "coio_C")


#### Process sample and gene level metadata ####
#Construct a design matrix from sample ids
design_matrix = constructDesignMatrix_SL1344(colnames(data[3:ncol(data)])) %>%
  dplyr::filter(!(donor == "fpdj" & replicate == 2)) #Remove the second replicate of fpdj

#Load sample metadata from disk
line_data = readRDS("annotations/compiled_line_metadata.rds") %>%
  #Keep only the lines that are also in the design matrix
  dplyr::semi_join(design_matrix, by = c("donor", "replicate")) 

#Load gene metadata from disk
length_df = dplyr::select(data, gene_id, length) %>% tbl_df()
gene_metadata = readRDS("annotations/biomart_transcripts.filtered.rds") %>% 
  dplyr::select(gene_id, gene_name, chromosome_name, gene_start, gene_end, percentage_gc_content) %>%
  dplyr::rename(chr = chromosome_name, left = gene_start, right = gene_end) %>%
  unique() %>% tbl_df() %>%
  dplyr::left_join(length_df, by = "gene_id")

#### Normalize gene expression data ####
#Construct a count matrix
count_matrix = data[,design_matrix$sample_id]
rownames(count_matrix) = data$gene_id

#Normalize gene expression data CQN
counts_selected = count_matrix[gene_metadata$gene_id,]
expression_cqn = cqn(counts = counts_selected, x = gene_metadata$percentage_gc_content, 
                     lengths = gene_metadata$length, verbose = TRUE)
expression_norm = expression_cqn$y + expression_cqn$offset

#Import genotype data from the VCF file
vcf_file = vcfToMatrix("genotypes/selected_genotypes.GRCh38.sorted.vcf.gz", "GRCh38")
snpspos = vcf_file$snpspos
genotypes = vcf_file$genotypes


#### Combine all of the data together ####
eqtl_dataset = list(design = design_matrix, line_metadata = line_data, 
                    gene_metadata = gene_metadata, exprs_counts = counts_selected,
                    exprs_cqn = expression_norm, snpspos = vcf_file$snpspos,
                    genotypes = vcf_file$genotypes)
saveRDS(eqtl_dataset,"results/SL1344/matrixeQTL/matrixeQTL_dataset.rds")



