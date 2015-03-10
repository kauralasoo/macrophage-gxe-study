library("MatrixEQTL")
library("ggplot2")
load_all("macrophage-gxe-study/seqUtils/")


#Load matrixeQTL dataset preapred previously by matrixeQTL_prepare.R script
eqtl_dataset = readRDS("results/SL1344/matrixeQTL/matrixeQTL_dataset.rds")

### Perform cis-eQTL mapping in the unstimulated condition ###

#Filter genes by expression
design_filtered = dplyr::filter(eqtl_dataset$design, condition == "D")
expression_filtered = eqtl_dataset$exprs_cqn[,design_filtered$sample_id]
expressed_genes = names(which(rowMeans(expression_filtered) > 2)) #Set conservative threshold to expression level
expression_sel = expression_filtered[expressed_genes,]

#Set up the genepos data.frame
genepos = dplyr::filter(eqtl_dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::select(gene_id, chr, left, right) %>% 
  dplyr::rename(geneid = gene_id) %>%
  as.data.frame()

#Perform PCA on the expression data
pca = performPCA(expression_sel, design_filtered)
ggplot(pca$pca_matrix, aes(x = PC1, y = PC2, label = sample_id)) + geom_point() + geom_text()
ggplot(pca$pca_matrix, aes(x = PC3, y = PC4, label = sample_id)) + geom_point() + geom_text()
ggplot(pca$pca_matrix, aes(x = PC5, y = PC6, label = sample_id)) + geom_point() + geom_text()


#Reconstruct expression data from first two PCs
pca_explained = (pca$pca_object$x[,1:4] %*% t(pca$pca_object$rotation[,1:4])) %>%
  scale(center = -1 * pca$pca_object$center, scale = FALSE) %>% t()
expression_sel = expression_sel - pca_explained

#Construct a SlicedData object of the expression data
expression_sliced = SlicedData$new()
expression_sliced$CreateFromMatrix(expression_sel)
expression_sliced$ResliceCombined(sliceSize = 2000)

#Match samples to genotyps
sample_meta = dplyr::left_join(design_filtered, eqtl_dataset$line_metadata, by = c("donor", "replicate"))
sample_genotype_match = dplyr::select(sample_meta, sample_id, genotype_id)
genotypes = eqtl_dataset$genotypes[,sample_genotype_match$genotype_id]
colnames(genotypes) = sample_genotype_match$sample_id

#Create a SlicedData obejct for the genotypes
snps = SlicedData$new()
snps$CreateFromMatrix(genotypes)
snps$ResliceCombined()

#Initialize matrixEQTL
output_file_name.cis = "results/SL1344/matrixeQTL/temp/matrixeQTL_condition_D.txt"

#RUN
me = Matrix_eQTL_main(
  snps = snps,
  gene = expression_sliced,
  cvrt = SlicedData$new(),
  output_file_name = "",
  pvOutputThreshold = 0,  
  output_file_name.cis = output_file_name.cis,
  pvOutputThreshold.cis = 1e-2,
  snpspos = eqtl_dataset$snpspos,
  genepos = genepos,
  cisDist = 5e5,
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name.cis);

cis_eqtls = dplyr::filter(me$cis$eqtls, FDR < 0.1)
length(unique(cis_eqtls$gene))
write.table(me$cis$eqtls, output_file_name.cis, sep ="\t", quote = FALSE, row.names = FALSE)






