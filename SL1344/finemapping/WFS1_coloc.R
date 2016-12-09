library("devtools")
library("plyr")
library("dplyr")
library("ggplot2")
library("purrr")
library("GenomicFeatures")
load_all("../seqUtils/")
load_all("../wiggleplotr")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)
rna_meta_str1 = wiggleplotrConstructMetadata(combined_expression_data$counts, combined_expression_data$sample_metadata, 
                                           "/Volumes/Ajamasin/bigwig/RNA/", bigWig_suffix = ".str1.bw") %>%
  dplyr::mutate(scaling_factor = 1)
rna_meta_str2 = wiggleplotrConstructMetadata(combined_expression_data$counts, combined_expression_data$sample_metadata, 
                                             "/Volumes/Ajamasin/bigwig/RNA/", bigWig_suffix = ".str2.bw") %>%
  dplyr::mutate(scaling_factor = 1)


#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, 
                                                  levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
atac_meta_df = wiggleplotrConstructMetadata(atac_list$counts, atac_list$sample_metadata, "/Volumes/Ajamasin/bigwig/ATAC/")


#### Import QTL mapping p-values ####
atac_rasqual_pvals = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
atac_fastqtl_pvals = readRDS("results/ATAC/QTLs/fastqtl_min_pvalues.rds")

rna_rasqual_pvals = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rna_fastqtl_pvals = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")


#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")


#### eQTL boxplot ####
#Make eQTL boxplot
WFS1_data = constructQtlPlotDataFrame("ENSG00000109501", "rs881796", combined_expression_data$cqn, vcf_file$genotypes, 
                                      combined_expression_data$sample_metadata, combined_expression_data$gene_metadata) %>%
  dplyr::filter(condition_name %in% c("naive","IFNg"))
WFS1_plot = plotQtlRow(WFS1_data)
ggsave("figures/main_figures/WFS1_RNA_boxplot.pdf", plot = WFS1_plot, width = 6, height = 4)


#Fetch Rasqual and fastqtl pvalues
#Fetch RASQUAL results for the ATAC peak
peak1_region = constructGeneRanges(data_frame(gene_id = "ATAC_peak_192294"), atac_list$gene_metadata, 1e5)
atac_rasqual_pvalues = purrr::map_df(qtlResults()$atac_rasqual, ~tabixFetchGenes(peak1_region, .)[[1]], .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
ggplot(dplyr::filter(atac_rasqual_pvalues, condition_name == "naive"), aes(x = pos, y = -log(p_nominal,10))) + geom_point()

peak1_region = constructGeneRanges(data_frame(gene_id = "ENSG00000109501"), combined_expression_data$gene_metadata, 1e5)
rna_rasqual_pvalues = purrr::map_df(qtlResults()$rna_rasqual, ~tabixFetchGenes(peak1_region, .)[[1]], .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
ggplot(dplyr::filter(rna_rasqual_pvalues, condition_name == "IFNg"), aes(x = pos, y = -log(p_nominal,10))) + geom_point()

#There seems to be another linked eQTL in the naive state (R2 = 0.48)
cor(vcf_file$genotypes["rs881796",],vcf_file$genotypes["rs33993436",] )



###### ATAC-seq coverage #####
#Construct metadata df for wiggleplotr
atac_track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs4689398", vcf_file$genotypes, -1) %>%
  dplyr::filter(track_id %in% c("naive", "IFNg"))

#Promoter region
region_coords = c(6269000, 6270500)
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 4, atac_list$gene_metadata)

#Make a coverage plot of the ATAC data
wfs1_promoter_atac = plotCoverage(
  exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
                              transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
                              connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
                              region_coords = region_coords, return_subplots_list = FALSE)
ggsave("figures/main_figures/WFS1_promoter_atac.pdf", plot = wfs1_promoter_atac, width = 3, height = 4)

#Enhancer region
region_coords = c(6285000, 6295000)
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 4, atac_list$gene_metadata)

#Make a coverage plot of the ATAC data
wfs1_enhancer_atac = plotCoverage(
  exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
  transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), 
  connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.7,0.3), 
  region_coords = region_coords, return_subplots_list = FALSE)
ggsave("figures/main_figures/WFS1_enhancer_atac.pdf", plot = wfs1_enhancer_atac, width = 4, height = 4)






#Plot gene structure
#Make a read coverage plot for RNA-Seq
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_79/TranscriptDb_GRCh38_79.db")
tx_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_79/Homo_sapiens.GRCh38.79.transcript_data.rds") %>%
  dplyr::rename(transcript_id = ensembl_transcript_id,
                gene_id = ensembl_gene_id,
                gene_name = external_gene_name)
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

#Filter transcripts for PTK2B
wfs1_tx = dplyr::filter(tx_metadata, gene_name == "WFS1", 
                         transcript_biotype == "protein_coding", transcript_status == "KNOWN")[2,]
region_coords = c(6269000,6303300)
region_coords = c(6250000,6320000)

tx_plot = plotTranscripts(exons[wfs1_tx$transcript_id], cdss[wfs1_tx$transcript_id], tx_metadata, rescale_introns = FALSE, 
                region_coords = region_coords)

#Fetch all peak in the region
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 4, atac_list$gene_metadata)
peak_plot = plotTranscripts(peak_annot$peak_list, peak_annot$peak_list, peak_annot$peak_annot, rescale_introns = FALSE, 
                               region_coords = region_coords, connect_exons = FALSE, label_type = "peak") + dataTrackTheme()


#Make a manhattan plot
atac = dplyr::filter(atac_rasqual_pvalues, condition_name %in% c("naive","IFNg")) %>% dplyr::mutate(track = paste(condition_name, "caQTL"))
rna = dplyr::filter(rna_rasqual_pvalues, condition_name %in% c("naive","IFNg")) %>% dplyr::mutate(track = paste(condition_name, "eQTL"))
pvals = bind_rows(atac, rna) %>% dplyr::arrange(p_nominal) %>% addR2FromLead(vcf_file$genotypes) %>%
  dplyr::mutate(track = factor(track, levels = c("naive caQTL", "IFNg caQTL", "naive eQTL", "IFNg eQTL")))



#Make a manhattan plot
rasqual_manhattan = ggplot(pvals, aes(x = pos, y = -log(p_nominal, 10), colour = R2)) + 
  geom_point() + 
  ylab(expression(paste("-",log[10], " p-value"))) +
  facet_grid(track~., scales = "free_y") +
  scale_x_continuous(expand = c(0,0), limits = region_coords) +
  theme_light() +
  dataTrackTheme()

joint_plot = cowplot::plot_grid(rasqual_manhattan, peak_plot, tx_plot, align = "v", ncol = 1, rel_heights = c(6,1,1.5))
ggsave("figures/main_figures/WFS1_manhattan.pdf", plot = joint_plot, width = 6, height = 8)


#Make read coverage plot of the gene
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2, "rs4689398", vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data, track_id %in% c("naive", "IFNg"))


WFS1_read_coverage = plotCoverage(exons = exons[wfs1_tx$transcript_id], cdss = cdss[wfs1_tx$transcript_id], track_data = filtered_tracks, rescale_introns = TRUE, 
             transcript_annotations = tx_metadata, fill_palette = getGenotypePalette(), 
             connect_exons = TRUE, label_type = "transcript", plot_fraction = 0.2, heights = c(0.7,0.3), 
             region_coords = region_coords, return_subplots_list = FALSE)
ggsave("figures/main_figures/WFS1_read_coverage.pdf", plot = WFS1_read_coverage, width = 6, height = 4)




#Import meQTL data from Pancreatic islets
me_variant
cg00701064
excel_colnames = c("snp_id","probe_id","tstat","pval","chr","pos", "snp_status","minor_allele", "MAF", 
                   "probe_pos","gene_names", "probe_desc")
me_summaries = readxl::read_excel("databases/methylation/journal.pgen.1004735.s029.XLSX", sheet = 1)
colnames(me_summaries) = excel_colnames
a = dplyr::filter(me_summaries, probe_id == "cg00701064")







#Check if the causal SNPs disrupt any known TF motifs
#Import motif matches
motif_metadata = readRDS("results/ATAC/cisBP/cisBP_motif_metadata.rds") %>%
  dplyr::transmute(motif_id = Motif_ID, tf_name = TF_Name, tf_count = TF_count)
motif_disruptions = importMotifDisruptions("results/ATAC/motif_analysis/motif_disruption.txt") %>%
  dplyr::left_join(motif_metadata, by = "motif_id")

#Filter by SNP ID
motif_hits = dplyr::filter(motif_disruptions, snp_id %in% c("rs7594476")) %>%
  dplyr::filter(max_rel_score > 0.8) %>% dplyr::arrange(-abs(rel_diff))


