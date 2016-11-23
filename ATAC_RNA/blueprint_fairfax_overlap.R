library("rhdf5")
library("dplyr")
library("tidyr")

h5f = H5Fopen("databases/Fairfax/hdf5/300_036.hdf5")
H5Fclose(h5f)

importBlueprintSummary <- function(summary_path){
  colnames = c("snp_string", "snp_id", "peak_id","p_nominal","beta", "p_bonferroni", "FDR", "alt_AF","std_error")
  data = readr::read_delim(summary_path, delim = "\t", col_names = colnames, col_types = "cccdddddd") %>%
    dplyr::filter(p_bonferroni <= 1) %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, "fdr")) %>%
    tidyr::separate(snp_string, c("chr","pos_string"), sep = ":") %>%
    tidyr::separate(pos_string, c("pos","ref","alt"), sep = "_") %>%
    dplyr::mutate(pos = as.integer(pos))
  return(data)
}

#Import BLUEPRINT min p-values
k4me1 = importBlueprintSummary("databases/BLUEPRINT/mono_K4ME1_min_pvalues.txt.gz")
k27ac = importBlueprintSummary("databases/BLUEPRINT/mono_K27AC_min_pvalues.txt.gz")

#Construct metadata for chromatin qtls
k4me1_meta = dplyr::select(k4me1, peak_id) %>% 
  tidyr::separate(peak_id, c("chr","start","end"), sep = ":", remove = FALSE, convert = TRUE)
k27ac_meta = dplyr::select(k27ac, peak_id) %>% 
  tidyr::separate(peak_id, c("chr","start","end"), sep = ":", remove = FALSE, convert = TRUE)

