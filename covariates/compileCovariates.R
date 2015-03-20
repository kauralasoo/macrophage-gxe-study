library("dplyr")
library("ggplot2")
library("tidyr")

#Import different files
dat = tbl_df(read.csv("macrophage-gxe-study/data/sample_lists/line_metadata_020315.csv", stringsAsFactors = FALSE, na.strings = ""))
flow_purity = readRDS("results/covariates/flow_cytometry_purity.rds")
rna_concentrations = readRDS("results/covariates/rna_concentrations.rds")

#Load genotype sample names
genotypes = read.table("genotypes/genotype_sample_names.txt", stringsAsFactors = FALSE)
colnames(genotypes) = c("genotype_id")
genotypes_names = tidyr::separate(genotypes, genotype_id, into = c("batch_id", "line_id"), sep = "-", remove = FALSE) %>% 
  dplyr::select(genotype_id, line_id)

#Convert all date fields to date
line_metadata = dat %>%
  tidyr::separate(line_id, into = c("donor","clone"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(ips_received = as.Date(ips_received, "%d/%m/%Y"), 
                ips_started = as.Date(ips_started, "%d/%m/%Y"),
                EB_formation = as.Date(EB_formation, "%d/%m/%Y"),
                diff_start = as.Date(diff_start, "%d/%m/%Y"),
                MF_harvest = as.Date(MF_harvest, "%d/%m/%Y"),
                salmonella = as.Date(salmonella, "%d/%m/%Y"),
                flow_date = as.Date(flow_date, "%d/%m/%Y"),
                terminated = as.Date(terminated, "%d/%m/%Y"),
                extraction_date = as.Date(extraction_date, "%d/%m/%Y"),
                rna_submit = as.Date(rna_submit, "%d/%m/%Y"))

#Add mean RNA concentration and flow purity to the data
mean_rna_concentrations = dplyr::select(rna_concentrations, donor, ng_ul_mean, replicate) %>% unique()
mean_flow_purity = flow_purity %>%
  dplyr::filter(!(donor == "gedo" & channel == "Pacific.Blue.A")) %>% #Remove an outlier measurement
  group_by(donor,flow_date) %>% 
  dplyr::summarize(mean_purity = mean(purity), max_purity = max(purity))

line_data = dplyr::left_join(line_metadata, mean_rna_concentrations, by = c("donor", "replicate")) %>%
  dplyr::left_join(mean_flow_purity, by = c("donor", "flow_date")) %>%
  dplyr::left_join(genotypes_names, by = "line_id")
saveRDS(line_data, "results/SL1344/sample_info/compiled_line_metadata.rds")

