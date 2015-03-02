library("flowCore")
library("flowViz")
library("tidyr")
library("dplyr")
source("macrophage-gxe-study/flow/flow_functions.R")

#Construct a flow metadata data.frame:
file_names = data.frame(file_name = list.files("flow/fcs/"), stringsAsFactors = FALSE)
metadata = tidyr::separate(file_names, file_name, c("name", "suffix"), sep = "\\.", remove = FALSE) %>% 
  tidyr::separate(name, into = c("donor", "date", "staining"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(sample = paste(donor, date, sep = "_")) %>% 
  dplyr::select(name, sample, donor, date, staining, file_name)

#Read all FCS files into a list
flow_files_list = as.list(paste("flow/fcs/", metadata$file_name, sep =""))
names(flow_files_list) = metadata$name
fcs_list = lapply(flow_files_list, read.FCS, alter.names = TRUE)

#Rename columns in FCS files
name_map = data.frame(old_names = c("X670.14..640..A","X450.50..405..A","X586.15..561..A"), 
                      new_names = c("APC.A","Pacific.Blue.A", "PE.A"), stringsAsFactors = FALSE)
fcs_list_renamed = lapply(fcs_list, changeFlowFrameNames, name_map)

#Keep only selected columns and merge into flowSet
selected_columns = c("FSC.A","SSC.A","APC.A","Pacific.Blue.A","PE.A","Time")
fcs_list_selected = lapply(fcs_list_renamed, function(flowframe, sel_cols){return(flowframe[,sel_cols])}, selected_columns)
flow_set = as(fcs_list_selected, "flowSet")
flow_set = addMetadataToFlowset(flow_set, metadata)

#Save the flow set to disk
saveRDS(flow_set, "flow/combined_flowSet.rds")



