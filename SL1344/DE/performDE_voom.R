library("limma")
library("edgeR")

#Load expression data from disk
expression_list = readRDS("results/SL1344/combined_expression_data.rds")

#Remove genes with very low counts
#nonzero_counts = expression_list$exprs_counts[apply(expression_list$exprs_counts, 1, mean) > 5,]

#Apply TMM normalisation
dge = DGEList(counts = expression_list$exprs_counts)
dge = calcNormFactors(dge, method = "TMM")

#Apply voom transformation
design_matrix = model.matrix(~SL1344 + IFNg + SL1344:IFNg, data = expression_list$design)
v <- voom(dge,design_matrix,plot=TRUE)
fit = lmFit(v, design_matrix)
fit <- eBayes(fit)
a = topTable(fit,coef=2)

inter = topTable(fit,coef=4, n = Inf)
ifng = topTable(fit,coef=3, n = Inf)
sl = topTable(fit,coef=2, n = Inf)

sl["ENSG00000175505",]
ifng["ENSG00000175505",]
inter["ENSG00000175505",]


