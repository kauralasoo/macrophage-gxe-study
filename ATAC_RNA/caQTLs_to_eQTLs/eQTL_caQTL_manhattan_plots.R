library("devtools")
library("plyr")
library("dplyr")
library("purrr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")
load_all("../wiggleplotr/")

#Functions
makeColocalisationPlot <- function(gene_id, peak_id, gene_metadata, peak_metadata, conditions){
  
  #Import ATAC-seq p-values from the region
  peak_pvalues = purrr::map_df(qtlResults()$atac_rasqual, ~tabixFetchGenesQuick(c(peak_id),
                                                                                tabix_file = ., 
                                                                                gene_metadata = peak_metadata, cis_window = 1e5)[[1]],
                               .id = "condition_name") %>%
    dplyr::filter(condition_name %in% conditions) %>%
    dplyr::mutate(condition_name = paste0("ATAC_",condition_name))
  peak_range = range(peak_pvalues$pos)
  
  #Import gene p-values from the region
  gene_pvalues = purrr::map_df(qtlResults()$rna_rasqual, ~tabixFetchGenesQuick(c(gene_id),
                                                                               tabix_file = ., 
                                                                               gene_metadata = gene_metadata, cis_window = 5e5)[[1]],
                               .id = "condition_name") %>%
    dplyr::filter(condition_name %in% conditions) %>%
    dplyr::mutate(condition_name = paste0("RNA_",condition_name)) %>%
    dplyr::filter(pos > peak_range[1], pos < peak_range[2])
  
  #Merge peak and gene pvalus
  factor_levels = c(paste0("ATAC_", conditions), paste0("RNA_", conditions))
  joint_pvalues = dplyr::bind_rows(peak_pvalues, gene_pvalues) %>%
    dplyr::arrange(p_nominal) %>%
    addR2FromLead(vcf_file$genotypes) %>%
    dplyr::mutate(track_id = factor(condition_name, levels = factor_levels))
  
  #Make Manhattan plots
  manhattan_plot = makeManhattanPlot(joint_pvalues, range(joint_pvalues$pos), color_R2 = TRUE, data_track = FALSE)
  return(manhattan_plot)
}

#Import RNA and ATAC data
rna_list = readRDS("results/SL1344/combined_expression_data_covariates.rds")
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import linked QTL pairs
pairs = readRDS("results/ATAC_RNA_overlaps/condition_specific_pairs.rds")

#makeColocalisationPlot("ENSG00000061918","ATAC_peak_204602",rna_list$gene_metadata, atac_list$gene_metadata, conditions = c("naive","IFNg"))

#IFNg pairs
ifng_pairs = dplyr::select(pairs$IFNg, gene_name, gene_id, peak_id) %>% unique()

plots = purrr::by_row(ifng_pairs, 
                      ~makeColocalisationPlot(.$gene_id, .$peak_id, rna_list$gene_metadata, atac_list$gene_metadata, conditions = c("naive","IFNg")), .to = "plot")
plot_list = setNames(plots$plot, paste(plots$gene_name, plots$peak_id, sep ="-"))
savePlotList(plot_list, "results/ATAC_RNA_overlaps/foreshadowing_plots/IFNg/", width = 6, height = 6)

#Salmonella pairs
sl1344_pairs = dplyr::select(pairs$SL1344, gene_name, gene_id, peak_id) %>% unique()

plots = purrr::by_row(sl1344_pairs, 
                      ~makeColocalisationPlot(.$gene_id, .$peak_id, rna_list$gene_metadata, atac_list$gene_metadata, conditions = c("naive","SL1344")), .to = "plot")
plot_list = setNames(plots$plot, paste(plots$gene_name, plots$peak_id, sep ="-"))
savePlotList(plot_list, "results/ATAC_RNA_overlaps/foreshadowing_plots/SL1344/", width = 6, height = 6)


