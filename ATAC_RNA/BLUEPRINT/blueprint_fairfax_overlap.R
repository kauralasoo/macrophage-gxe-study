library("rhdf5")
library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
load_all("../seqUtils/")

importBlueprintSummary <- function(summary_path){
  colnames = c("snp_string", "old_snp_id", "peak_id","p_nominal","beta", "p_bonferroni", "FDR", "alt_AF","std_error")
  data = readr::read_delim(summary_path, delim = "\t", col_names = colnames, col_types = "cccdddddd") %>%
    dplyr::filter(p_bonferroni <= 1) %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, "fdr")) %>%
    tidyr::separate(snp_string, c("chr","pos_string"), sep = ":") %>%
    tidyr::separate(pos_string, c("pos","ref","alt"), sep = "_") %>%
    dplyr::mutate(pos = as.integer(pos))
  return(data)
}

extractProbeSummaries <- function(selected_probe_id, selected_snp_ids, hdf5_file){
  
  #Read probe data from the hdf5 file
  probe_data = rhdf5::h5read(hdf5_file, selected_probe_id)
  
  #Extracte relevant data from the hdf5 object
  conditions = probe_data$lmm$context
  snp_ids = probe_data$snp_info$genotype_id
  beta_cols = paste0(conditions, "_beta")
  pvalue_cols = paste0(conditions, "_pvalue")
  nsnps = probe_data$lmm$nusnps
  
  #Extract betas
  beta_matrix = as.data.frame(probe_data$lmm$beta) %>% tbl_df()
  colnames(beta_matrix) = beta_cols
  beta_df = dplyr::mutate(beta_matrix, old_snp_id = snp_ids) %>% 
    dplyr::filter(old_snp_id %in% selected_snp_ids)
  
  #Extract p_values
  p_matrix = as.data.frame(probe_data$lmm$pv) %>% tbl_df()
  colnames(p_matrix) = pvalue_cols
  p_df = dplyr::mutate(p_matrix, old_snp_id = snp_ids) %>% 
    dplyr::filter(old_snp_id %in% selected_snp_ids)
  
  #Merge results
  result = dplyr::left_join(beta_df, p_df, by = "old_snp_id") %>% 
    dplyr::mutate(n_snps = nsnps, probe_id = selected_probe_id) %>%
    dplyr::select(probe_id, old_snp_id, n_snps, everything())
  return(result)
}

extractBatchSummaries <- function(batch_df, hdf5_folder){
  
  #Add second batch_id
  batch_df = dplyr::mutate(batch_df, probe_id2 = probe_id)
  
  #Open the hdf5 file
  file_name = batch_df$file_name[1]
  print(file_name)
  file_path = file.path(hdf5_folder, file_name)
  h5f = rhdf5::H5Fopen(file_path)
  
  #Split batch into probe slices
  probe_slices = dplyr::group_by(batch_df, probe_id2) %>% 
    purrr::by_slice(identity)
  probe_slices_list = probe_slices$.out
  
  #Extract summaries for each slices
  batch_summaries = purrr::map_df(probe_slices_list, ~extractProbeSummaries(.$probe_id[1], .$old_snp_id, h5f))
  
  #Close hdf5 file
  rhdf5::H5Fclose(h5f)
  
  return(batch_summaries)
}

#Import BLUEPRINT min p-values
k4me1 = importBlueprintSummary("databases/BLUEPRINT/mono_K4ME1_min_pvalues.txt.gz")
k27ac = importBlueprintSummary("databases/BLUEPRINT/mono_K27AC_min_pvalues.txt.gz")

#Import GRCh37 variant information
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

#Filter k4me1 and k27ac results for SNPs that have genotype data
new_snps = dplyr::select(GRCh37_variants, chr, pos, snp_id) %>% dplyr::filter(pos %in% k4me1$pos)
k4me1_filtered = dplyr::left_join(k4me1, new_snps, by = c("chr","pos")) %>% dplyr::filter(!is.na(snp_id))

new_snps = dplyr::select(GRCh37_variants, chr, pos, snp_id) %>% dplyr::filter(pos %in% k27ac$pos)
k27ac_filtered = dplyr::left_join(k27ac, new_snps, by = c("chr","pos")) %>% dplyr::filter(!is.na(snp_id))

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Construct metadata for chromatin qtls
k4me1_meta = dplyr::select(k4me1_filtered, peak_id) %>% 
  tidyr::separate(peak_id, c("chr","start","end"), sep = ":", remove = FALSE, convert = TRUE)
k27ac_meta = dplyr::select(k27ac_filtered, peak_id) %>% 
  tidyr::separate(peak_id, c("chr","start","end"), sep = ":", remove = FALSE, convert = TRUE)

#Import Fairfax lead variant results
fairfax_summary = readr::read_delim("databases/Fairfax/summstat_1Mb.csv", delim = ",")
gene_metadata = fairfax_summary[,34:43] %>% 
  tidyr::separate(file, c("folder", "file_name"), sep = "0000/") %>% 
  dplyr::transmute(gene_id = geneID, chr = gene_chrom, start = gene_start, end = gene_end, 
                   length = gene_len, strand = gene_strand, tss = gene_tss, n_snps = nusnps, 
                   probe_id = probeID, file_name)
gene_probe_map = dplyr::select(gene_metadata, gene_id, probe_id)

#Extract p-vals for each condition
naive = dplyr::transmute(fairfax_summary, probe_id = probeID, chr = gene_chrom, 
                        pos = CD14_pos, snp_id = CD14_rs, ref = CD14_ref, alt = CD14_alt, 
                        p_nominal = CD14_pv_raw, p_bonferroni = CD14_pv_bonf, 
                        beta = CD14_beta, tstat = CD14_tstat)
IFN24 = dplyr::transmute(fairfax_summary, probe_id = probeID, chr = gene_chrom, 
                        pos = IFN_pos, snp_id = IFN_rs, ref = IFN_ref, alt = IFN_alt, 
                        p_nominal = IFN_pv_raw, p_bonferroni = IFN_pv_bonf, 
                        beta = IFN_beta, tstat = IFN_tstat)
LPS2 = dplyr::transmute(fairfax_summary, probe_id = probeID, chr = gene_chrom, 
                       pos = LPS2_pos, snp_id = LPS2_rs, ref = LPS2_ref, alt = LPS2_alt, 
                       p_nominal = LPS2_pv_raw, p_bonferroni = LPS2_pv_bonf, 
                       beta = LPS2_beta, tstat = LPS2_tstat)
LPS24 = dplyr::transmute(fairfax_summary, probe_id = probeID, chr = gene_chrom, 
                        pos = LPS24_pos, snp_id = LPS24_rs, ref = LPS24_ref, alt = LPS24_alt, 
                        p_nominal = LPS24_pv_raw, p_bonferroni = LPS24_pv_bonf, 
                        beta = LPS24_beta, tstat = LPS24_tstat)
fairfax_conditions = list(naive = naive, IFN24 = IFN24, LPS2 = LPS2, LPS24 = LPS24)

#Add FDR to each fo the conditions
fairfax_fdr = purrr::map(fairfax_conditions, ~dplyr::mutate(., p_bonferroni = pmin(p_bonferroni, 1)) %>%
                           dplyr::mutate(p_fdr = p.adjust(p_bonferroni, "fdr")))

#Identify significant QTLs in each condition
fairfax_hits = purrr::map(fairfax_fdr, ~dplyr::filter(., p_fdr < 0.1)) %>%
  purrr::map_df(identity, .id = "condition_name") %>%
  dplyr::left_join(gene_probe_map, by = "probe_id") %>%
  dplyr::group_by(condition_name, gene_id) %>%
  dplyr::arrange(condition_name, gene_id, p_nominal) %>%
  dplyr::filter(row_number() == 1) %>% 
  ungroup() %>%
  dplyr::select(gene_id, probe_id, snp_id, p_nominal) %>% 
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::select(gene_id, probe_id, snp_id) %>%
  dplyr::rename(old_snp_id = snp_id) %>%
  unique()

#Identify all lead variants
fairfax_lead_variants = purrr::map_df(fairfax_fdr, identity) %>% 
  dplyr::select(snp_id, chr, pos) %>% unique() %>% 
  dplyr::rename(old_snp_id = snp_id) %>%
  dplyr::mutate(chr = as.character(chr))

#Match them to macrophage genotype data
new_snp_ids = dplyr::select(GRCh37_variants, chr, pos, snp_id) %>% 
  dplyr::filter(pos %in% fairfax_lead_variants$pos)
matched_snp_ids = dplyr::left_join(fairfax_lead_variants, new_snp_ids, by = c("chr","pos")) %>% 
  dplyr::filter(!is.na(snp_id))

#Identify unique gene-snp pairs (based on R2)
fairfax_pairs = dplyr::left_join(fairfax_hits, matched_snp_ids, by = "old_snp_id") %>%
  dplyr::filter(!is.na(snp_id)) %>% 
  dplyr::select(gene_id, snp_id) %>% unique()
fairfax_filtered_pairs = filterHitsR2(fairfax_pairs, vcf_file$genotypes, R2_thresh = 0.8)

#Go back to probes
fairfax_filtered_hits = dplyr::left_join(fairfax_filtered_pairs, 
                                         dplyr::select(fairfax_hits, gene_id, probe_id) %>% unique(), by = "gene_id") %>%
  dplyr::left_join(matched_snp_ids, by = "snp_id")

#Extract probe-level summary stats for each gene-snp pair in each condition
hdf5_files = list.files("databases/Fairfax/hdf5")
fairfax_files = dplyr::left_join(fairfax_filtered_hits, 
                dplyr::select(gene_metadata, probe_id, file_name), by = "probe_id") %>%
  dplyr::arrange(file_name) %>% 
  dplyr::filter(file_name %in% hdf5_files)

#Iterate over HDF5 files
file_list = (dplyr::mutate(fairfax_files, file_name2 = file_name) %>% 
               dplyr::group_by(file_name2) %>% 
               purrr::by_slice(identity))$.out

#Extract summary stats for all variants
summaries = purrr::map_df(file_list, ~extractBatchSummaries(., hdf5_folder = "databases/Fairfax/hdf5/"))
saveRDS(summaries, "databases/Fairfax/snp_level_summaries.rds")
summaries = readRDS("databases/Fairfax/snp_level_summaries.rds")

#Add metadata
summaries_df = dplyr::select(fairfax_files, -file_name) %>%
  dplyr::left_join(summaries, by = c("probe_id","old_snp_id"))

#Calculate bonferroni pvalues
summaries_df_bonferroni = dplyr::mutate(summaries_df, CD14_pvalue = pmin(CD14_pvalue*n_snps,1),
              IFN_pvalue = pmin(IFN_pvalue*n_snps,1),
              LPS2_pvalue = pmin(LPS2_pvalue*n_snps,1),
              LPS24_pvalue = pmin(LPS24_pvalue*n_snps,1))

#Extract condition-specific eQTLs
IFN_eqtls = dplyr::filter(summaries_df_bonferroni, 
                          IFN_pvalue < 0.01, CD14_pvalue > 0.1, 
                          abs(CD14_beta) < 0.5, abs(IFN_beta) > 0.7)
LPS2_eqtls = dplyr::filter(summaries_df_bonferroni, 
                          LPS2_pvalue < 0.01, CD14_pvalue > 0.1, 
                          abs(CD14_beta) < 0.5, abs(LPS2_beta) > 0.7)
LPS24_eqtls = dplyr::filter(summaries_df_bonferroni, 
                           LPS24_pvalue < 0.01, CD14_pvalue > 0.1, 
                           abs(CD14_beta) < 0.5, abs(LPS24_beta) > 0.7)
CD14_eqtls = dplyr::filter(summaries_df_bonferroni, CD14_pvalue < 0.01, abs(CD14_beta) > 0.7)

other_eqtls =  dplyr::filter(summaries_df_bonferroni, 
                          CD14_pvalue > 0.1, 
                          abs(CD14_beta) < 0.5, 
                          pmax(abs(IFN_beta), abs(LPS2_beta), abs(LPS24_beta)) > 0.7,
                          min(IFN_pvalue, LPS2_pvalue, LPS24_pvalue) < 0.01)

CD14_not_eqtls = dplyr::filter(summaries_df_bonferroni, CD14_pvalue > 0.1, abs(CD14_beta) < 0.5)


#FInd overlaps between K4ME1 QTLs and eQTLs
k4me1_qtls = dplyr::filter(k4me1_filtered, p_fdr < 0.1)
k4me1_overlaps = findGWASOverlaps(dplyr::transmute(summaries_df_bonferroni, gene_id = probe_id, snp_id), 
                                     dplyr::select(k4me1_qtls, peak_id, snp_id, chr, pos), 
                                  vcf_file, max_distance = 5e5, min_r2 = 0.8)

k27ac_qtls = dplyr::filter(k27ac_filtered, p_fdr < 0.1)
k27ac_overlaps = findGWASOverlaps(dplyr::transmute(summaries_df_bonferroni, gene_id = probe_id, snp_id), 
                                  dplyr::select(k27ac_qtls, peak_id, snp_id, chr, pos), 
                                  vcf_file, max_distance = 5e5, min_r2 = 0.8)


k4me1_olaps_filtered = dplyr::rename(k4me1_overlaps, probe_id = gene_id) %>% 
  dplyr::filter(R2 > 0.8) %>%
  dplyr::select(probe_id, snp_id) %>%
  unique()
k27ac_olaps_filtered = dplyr::rename(k27ac_overlaps, probe_id = gene_id) %>% 
  dplyr::filter(R2 > 0.8) %>%
  dplyr::select(probe_id, snp_id) %>%
  unique()

dplyr::semi_join(k4me1_olaps_filtered, IFN_eqtls, by = c("probe_id", "snp_id"))
dplyr::semi_join(k27ac_olaps_filtered, IFN_eqtls, by = c("probe_id", "snp_id"))

dplyr::semi_join(k4me1_olaps_filtered, CD14_eqtls, by = c("probe_id", "snp_id"))
dplyr::semi_join(k27ac_olaps_filtered, CD14_eqtls, by = c("probe_id", "snp_id"))

dplyr::semi_join(k4me1_olaps_filtered, other_eqtls, by = c("probe_id", "snp_id"))
dplyr::semi_join(k27ac_olaps_filtered, other_eqtls, by = c("probe_id", "snp_id"))


dplyr::semi_join(k4me1_olaps_filtered, LPS2_eqtls, by = c("probe_id", "snp_id"))
dplyr::semi_join(k27ac_olaps_filtered, LPS2_eqtls, by = c("probe_id", "snp_id"))

dplyr::semi_join(k4me1_olaps_filtered, LPS24_eqtls, by = c("probe_id", "snp_id"))
dplyr::semi_join(k27ac_olaps_filtered, LPS24_eqtls, by = c("probe_id", "snp_id"))

