#Import coordinates
snp_coords = readr::read_delim("genotypes/SL1344/imputed_20151005/imputed.86_samples.snp_coords.txt", 
                               delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))
union_exon_coords = read.table("annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.filtered_genes.txt", 
                               stringsAsFactors = FALSE, header = TRUE) %>% dplyr::rename(chr = chromosome_name)
filtered_coords = dplyr::semi_join(union_exon_coords, rna_conditions_renamed$naive$gene_metadata, by = "gene_id")
exon_df = countSnpsOverlapingExons(filtered_coords, snp_coords, cis_window = 500000) %>% dplyr::arrange(chromosome_name, range_start)


rasqual_input_folder = "results/SL1344/rasqual/input/"


### Identify which genes failed and construct new batches for them
#Naive
completed_genes = readr::read_delim("results/SL1344/rasqual/output/naive_500kb/naive.completed_genes.txt", col_names = "gene_id", delim = "\t")
failed_genes = dplyr::anti_join(exon_df, completed_genes)
failed_batches = rasqualConstructGeneBatches(failed_genes, 1, "batch_5")
write.table(failed_batches, file.path(rasqual_input_folder, "naive.failed_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#IFNg
completed_genes = readr::read_delim("results/SL1344/rasqual/output/IFNg_500kb/IFNg.completed_genes.txt", col_names = "gene_id", delim = "\t")
failed_genes = dplyr::anti_join(exon_df, completed_genes, by = "gene_id")
failed_batches = rasqualConstructGeneBatches(failed_genes, 1, "batch_5")
write.table(failed_batches, file.path(rasqual_input_folder, "IFNg.failed_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#SL1344
completed_genes = readr::read_delim("results/SL1344/rasqual/output/SL1344_500kb/SL1344.completed_genes.txt", col_names = "gene_id", delim = "\t")
failed_genes = dplyr::anti_join(exon_df, completed_genes, by = "gene_id")
failed_batches = rasqualConstructGeneBatches(failed_genes, 1, "batch_5")
write.table(failed_batches, file.path(rasqual_input_folder, "SL1344.failed_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#SL1344
completed_genes = readr::read_delim("results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344.completed_genes.txt", col_names = "gene_id", delim = "\t")
failed_genes = dplyr::anti_join(exon_df, completed_genes, by = "gene_id")
failed_batches = rasqualConstructGeneBatches(failed_genes, 1, "batch_5")
write.table(failed_batches, file.path(rasqual_input_folder, "IFNg_SL1344.failed_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


