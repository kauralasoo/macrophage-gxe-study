library("DESeq2")
library("dplyr")

#Combine counts
sample_names = read.table("fastq/acLDL_samples.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,2]

loadCounts <- function(sample_dir, sample_names, counts_suffix = ".counts.txt" ){
  matrix = c()
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    print(path)
    table = read.table(path, header = TRUE)
    print(head(table))
    if (i == 1){
      matrix = table[,c(1,6,7)]
    }
    else{
      matrix = cbind(matrix, table[,7])
    }
  }
  colnames(matrix) = c("gene_id", "length", sample_names)
  return(matrix)
}

data = loadCounts("STAR/", sample_names)

#Separate count matrix from data
count_matrix = data
rownames(count_matrix) = data$gene_id
count_matrix = count_matrix[,3:ncol(count_matrix)]

#Construct a design matrix
design = data.frame(sample_id = sample_names, treatment = c("ctrl", "acLDL", "ctrl", "acLDL", "ctrl", "acLDL"))
rownames(design) = design$sample_id

dds = DESeqDataSetFromMatrix(count_matrix, design, ~treatment)
dds = DESeq(dds)
res = results(dds)
de_genes = dplyr::filter(res, padj < 0.05, abs(log2FoldChange) > 0.58) %>% dplyr::arrange(desc(log2FoldChange))
