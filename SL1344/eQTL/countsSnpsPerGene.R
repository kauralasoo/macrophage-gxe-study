library("GenomicRanges")
library("plyr")
library("dplyr")

#Count the number of cis SNPs and feature SNPs for each gene
snp_coords = read.table("genotypes/SL1344/snp_coords.txt", stringsAsFactors = FALSE)
snp_granges = GRanges(seqnames = snp_coords[,1], ranges = IRanges(start = snp_coords[,2], end = snp_coords[,2]))

#Import exon coordinates
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
union_exon_coords = read.table("annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.txt", stringsAsFactors = FALSE)
colnames(union_exon_coords) = c("gene_id", "exon_starts", "exon_ends")
union_exon_coords = union_exon_coords[union_exon_coords$gene_id %in% expression_dataset$gene_metadata$gene_id,]
gene_meta = dplyr::select(expression_dataset$gene_metadata, gene_id, chromosome_name, strand)
union_exon_coords = dplyr::left_join(union_exon_coords, gene_meta, by = "gene_id") %>%
  dplyr::arrange(chromosome_name)

#Create a df of exon coords
exon_df = ldply(dlply(union_exon_coords, .(gene_id), 
    function(x){
        data.frame(gene_id = x$gene_id, 
                        chromosome_name = x$chromosome_name,
                        strand = x$strand,
                        exon_start = as.numeric(unlist(strsplit(x$exon_starts,","))),
                        exon_end = as.numeric(unlist(strsplit(x$exon_ends,",")))
        )
      }
    )
  )
exon_df$gene_id = as.character(exon_df$gene_id)
exon_granges = GRanges(seqnames = exon_df$chromosome_name, ranges = IRanges(start = exon_df$exon_start, end = exon_df$exon_end), strand = exon_df$strand)
mcols(exon_granges) = data.frame(gene_id = exon_df$gene_id, stringsAsFactors = FALSE)

#Count the number of feature SNPS per gene
olaps = findOverlaps(exon_granges, snp_granges) %>% as.data.frame()
olaps = dplyr::mutate(olaps, gene_id = exon_granges$gene_id[olaps$queryHits])
feature_snp_count = dplyr::group_by(olaps, gene_id) %>% dplyr::summarise(feature_snp_count = length(subjectHits)) %>%
  dplyr::mutate(gene_id = as.character(gene_id))

#Count the number of cis SNPs per gene
gene_window = dplyr::group_by(exon_df, gene_id) %>% summarise(start = min(exon_start), end = max(exon_end))

#500 kb around TSS
TSS_window = dplyr::group_by(exon_df, gene_id) %>% 
  summarise(chromosome_name = chromosome_name[1], start = ifelse(max(strand) == 1, min(exon_start), max(exon_end)))
TSS_granges = GRanges(seqnames = TSS_window$chromosome_name, 
                      IRanges(start = pmax(TSS_window$start - 500000, 0), end = TSS_window$start + 500000))
mcols(TSS_granges) = data.frame(gene_id = TSS_window$gene_id)

TSS_olaps = findOverlaps(TSS_granges, snp_granges) %>% as.data.frame()
TSS_olaps = dplyr::mutate(TSS_olaps, gene_id = as.character(TSS_granges$gene_id)[TSS_olaps$queryHits])
TSS_cis_snps_count = dplyr::group_by(TSS_olaps, gene_id) %>% dplyr::summarise(cis_snp_count = length(subjectHits))

TSS_df = as.data.frame(TSS_granges) %>% tbl_df() %>% 
  dplyr::transmute(gene_id = as.character(gene_id), chr = seqnames, range_start = start, range_end = end) %>% 
  dplyr::left_join(feature_snp_count, by = "gene_id") %>% 
  dplyr::left_join(TSS_cis_snps_count, by = "gene_id")
TSS_df[is.na(TSS_df)] = 0

write.table(TSS_df, "rasqual/input/snp_counts.txt", row.names = FALSE, sep ="\t", quote = FALSE)

