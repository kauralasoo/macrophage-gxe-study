library("readr")
library("dplyr")
library("tidyr")
library("limma")
library("purrr")
load_all("../seqUtils/")

#Import read counts
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
sample_ids = c("chrom", combined_expression_data$sample_metadata$sample_id)

#Import intron and cluster counts from leafcutter
intron_counts = readr::read_delim("results/SL1344/leafcutter/leafcutter.intron_counts.txt.gz", 
                                  delim = " ", col_names = TRUE)[,sample_ids]
cluster_counts = readr::read_delim("results/SL1344/leafcutter/leafcutter.cluster_counts.txt.gz", 
                                  delim = " ", col_names = TRUE)[,sample_ids]

#Convert to matrices
intron_matrix = dplyr::select(intron_counts, -chrom) %>% as.matrix()
rownames(intron_matrix) = intron_counts$chrom

cluster_matrix = dplyr::select(cluster_counts, -chrom) %>% as.matrix()
rownames(cluster_matrix) = cluster_counts$chrom

#Filter introns
filtered_introns = intron_matrix[rowSums(intron_matrix >= 1) >= (ncol(intron_counts) - 20),]
filtered_clusters = cluster_matrix[rownames(filtered_introns),]

#Extract intron coords from counts
intron_coords = dplyr::filter(intron_counts, chrom %in% rownames(filtered_introns)) %>%
  dplyr::select(chrom) %>%
  tidyr::separate(chrom, c("chr","start","end","cluster_id"), sep = ":", remove = FALSE) %>%
  dplyr::mutate(start = as.integer(start), end = as.integer(end)) %>%
  dplyr::rename(gene_id = chrom)

#Cluster metadata
cluster_meta = dplyr::group_by(intron_coords, cluster_id) %>%
  dplyr::summarise(chr = chr[1], start = min(start), end = max(end), cluster_size = length(cluster_id))

#Junction metadata
junction_meta = dplyr::select(intron_coords, gene_id, cluster_id) %>% 
  dplyr::left_join(cluster_meta, by = "cluster_id") %>%
  dplyr::mutate(strand = 1)

#Normalise intron exclusion proportions
intron_protortions = filtered_introns/filtered_clusters
intron_protortions[is.nan(intron_protortions)] = 0

#Standardise intron proportions
intron_prop_std = ExpressionSet(intron_protortions) %>% standardise() %>% exprs()
intron_prop_quantile = quantileNormaliseMatrix(intron_prop_std)

#Calculate covariates
sample_metadata = combined_expression_data$sample_metadata[,-grep("PEER*", colnames(combined_expression_data$sample_metadata))]

#Perform PCA on the intron proportions
a = performPCA(intron_prop_quantile, sample_metadata)
ggplot(a$pca_matrix, aes(x = PC1, y = PC2, color = rna_auto)) + geom_point()

#Make a list of prop data
prop_list = list(counts = filtered_introns,
               tpm = intron_protortions,
               cqn = intron_prop_quantile,
               norm_factors = combined_expression_data$norm_factors,
               sample_metadata = sample_metadata,
               gene_metadata = junction_meta)

#Extraxt conditions from the list
condition_list = list(naive = "naive",SL1344 = "SL1344",IFNg = "IFNg", IFNg_SL344 = "IFNg_SL1344")
prop_conditions = purrr::map(condition_list, ~extractConditionFromExpressionList(.,prop_list))

#Use PCA on nromalised proportion values to construct covariates for QTL mapping
prop_pca = purrr::map(prop_conditions, ~performPCA(.$cqn, .$sample_metadata, n_pcs = 10, column_prefix = "norm_"))
ggplot(prop_pca$naive$pca_matrix, aes(x = norm_PC1, y = norm_PC2, color = rna_auto)) + geom_point()

#Add PCs into the sample metadata df
covariates = purrr::map(prop_pca, ~.$pca_matrix) %>% 
  dplyr::bind_rows() %>%
  dplyr::mutate(sex_binary = ifelse(sex == "male",0,1))
new_sample_metadata = dplyr::select(prop_list$sample_metadata, sample_id) %>% dplyr::left_join(covariates, by = "sample_id")
prop_list$sample_metadata = new_sample_metadata

#Export data to disk
saveRDS(prop_list, "results/SL1344/combined_proportions.rds")
