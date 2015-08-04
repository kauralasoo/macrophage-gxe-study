library("plyr")
library("dplyr")
library("devtools")
library("rtracklayer")
load_all("macrophage-gxe-study/seqUtils/")

#Import exon start-end coordinates
exon_start_end = read.table("annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.txt", stringsAsFactors = FALSE)
colnames(exon_start_end) = c("gene_id", "exon_starts", "exon_ends")

#Import gene metadata
metadata = readRDS("annotations/Homo_sapiens.GRCh38.79.transcript_data.rds")
gene_metadata = dplyr::select(metadata, ensembl_gene_id, gene_biotype, chromosome_name, strand) %>% 
  unique() %>%
  dplyr::rename(gene_id = ensembl_gene_id, chr = chromosome_name)

valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")

#Compile gene data into a data frame
gene_data = dplyr::left_join(exon_start_end, gene_metadata, by = "gene_id") %>%
  dplyr::filter(gene_biotype %in% valid_gene_biotypes, chr %in% valid_chromosomes)

#Make a list of gene ids
gene_id_list = as.list(gene_data$gene_id)
names(gene_id_list) = gene_data$gene_id

#Make gff file of intron coorinates
intron_df_list = lapply(gene_id_list, constructIntronExonDf, gene_data, intron_gap = 10, type = "intron")
intron_df_list = intron_df_list[!unlist(lapply(intron_df_list, function(x){is.null(x)}))]#remove nulls
intron_df = ldply(intron_df_list) %>% dplyr::select(-.id) %>% dplyr::filter(end - start > 0)
saveRDS(intron_df, "annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.intron_df.rds")
intron_gr = dataFrameToGRanges(intron_df)
export.gff3(intron_gr, "annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.introns.gff3")

#Make gff file of exon coordinates
exon_df_list = lapply(gene_id_list, constructIntronExonDf, gene_data, intron_gap = 10, type = "exon")
intron_genes = unique(intron_df$gene_id)
exon_df = ldply(exon_df_list) %>% dplyr::select(-.id) %>% dplyr::filter(gene_id %in% intron_genes)
saveRDS(exon_df, "annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.exon_df.rds")
exon_gr = dataFrameToGRanges(exon_df)
export.gff3(exon_gr, "annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.exons.gff3")

