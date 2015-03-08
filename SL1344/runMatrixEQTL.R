library("dplyr")
library("devtools")
library("cqn")
load_all("macrophage-gxe-study/seqUtils/")

#Load raw counts from disk
data = read.table("results//SL1344//combined_counts.txt", stringsAsFactors = FALSE) #Read counts
#There seems to be a sample swap between coio_C and coio_D, fix that
indexes = c(which(colnames(data) == "coio_C"), which(colnames(data) == "coio_D"))
colnames(data)[indexes] = c("coio_D", "coio_C")

#Construct a design matrix from sample ids
design_matrix = constructDesignMatrix_SL1344(colnames(data[3:ncol(data)])) %>%
  dplyr::filter(!(donor == "fpdj" & replicate == 2)) #Remove the second replicate of fpdj

#Construct a count matrix
count_matrix = data[,design_matrix$sample_id]
rownames(count_matrix) = data$gene_id

#Load sample metadata from disk
line_data = readRDS("annotations/compiled_line_metadata.rds") %>%
  #Keep only the lines that are also in the design matrix
  dplyr::semi_join(design_matrix, by = c("donor", "replicate")) 

#Load gene metadata from disk
length_df = dplyr::select(data, gene_id, length) %>% tbl_df()
gene_metadata = readRDS("annotations/biomart_transcripts.filtered.rds") %>% 
  dplyr::select(gene_id, chromosome_name, gene_start, gene_end, percentage_gc_content) %>%
  dplyr::rename(chr = chromosome_name, left = gene_start, right = gene_end) %>%
  unique() %>% tbl_df()
gene_covariates = dplyr::left_join(gene_metadata, length_df, by = "gene_id")

#Normalize gene expression data CQN
counts_selected = count_matrix[gene_metadata$gene_id,]
expression_cqn = cqn(counts = counts_selected, x = gene_covariates$percentage_gc_content, 
                     lengths = gene_covariates$length, verbose = TRUE)
expression_norm = expression_cqn$y + expression_cqn$offset
saveRDS(expression_norm, "results/SL1344/cqn_normalized_expression.rds")
expression_norm = readRDS("results/SL1344/cqn_normalized_expression.rds")

### Perform cis-eQTL mapping in the unstimulated condition ###
#Filter genes by expression
design_filtered = dplyr::filter(design_matrix, condition == "A")
expression_filtered = expression_norm[,design_filtered$sample_id]
expressed_genes = names(which(rowMeans(expression_filtered) > 2))

#Construct a sample metadata matrix
sample_meta = dplyr::left_join(design_filtered, line_data, by = c("donor", "replicate"))

#Set up the gene expression matrix
expression_file_name = "results/SL1344/matrixEQTL/temp/expression_matrix.txt"
expression_matrix = expression_filtered[expressed_genes,]
write.table(expression_matrix, expression_file_name, sep ="\t", quote = FALSE)
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

#Set up the genepos data.frame
genepos = dplyr::filter(gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::select(gene_id, chr, left, right) %>% 
  dplyr::rename(geneid = gene_id) %>%
  as.data.frame()
snpspos = readRDS("genotypes/snp_positions.rds") %>% data.frame()

#Load genotype matrix from disk
genotype_file_name = "results/SL1344/matrixEQTL/temp/genotype_matrix.txt" 
genotypes = readRDS("genotypes/genotypes_for_matrixeQTL.rds")

#Match samples to genotypes
sample_genotype_match = dplyr::select(sample_meta, sample_id, genotype_id)
genotypes_reorder = genotypes[,sample_genotype_match$genotype_id]
colnames(genotypes_reorder) = sample_genotype_match$sample_id
write.table(genotypes_reorder, genotype_file_name, sep ="\t", quote = FALSE)

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(genotype_file_name)

#Initialize matrixEQTL
output_file_name.cis = "results/SL1344/matrixEQTL/temp/matrixeQTL_cis_run_output.txt"

#RUN
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = SlicedData$new(),
  output_file_name = "",
  pvOutputThreshold = 0,  
  output_file_name.cis = output_file_name.cis,
  pvOutputThreshold.cis = 1e-2,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = 1e6,
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name.cis);

cis_eqtls = dplyr::filter(me$cis$eqtls, FDR < 0.1)









