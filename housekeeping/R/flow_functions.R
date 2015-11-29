estimateSamplePurity <- function(sample_id, metadata, flowset, staining = "CD14+CD16+CD206", control = "isotype", channels = c("APC.A", "PE.A", "Pacific.Blue.A")){
  
  #Extraxt sample names from metadata
  stained_sample_name = dplyr::filter(filtered_metadata, sample == sample_id, staining == "CD14+CD16+CD206")$name
  control_sample_name = dplyr::filter(filtered_metadata, sample == sample_id, staining == "isotype")$name
  
  #Extract sample data
  stained = flowset[[stained_sample_name]] %>% exprs() %>% as.data.frame() %>% tbl_df()
  control = flowset[[control_sample_name]] %>% exprs() %>% as.data.frame() %>% tbl_df()
  combined = rbind(stained, control)
  
  results = c()
  #Iterate over all channels
  for (channel in channels){
    #Fit Mclust with either 2 or 3 components and choose the best number
    model2 = Mclust(combined[,channel], modelNames = "V", G = 2)
    model3 = Mclust(combined[,channel], modelNames = "V", G = 3)
    fraction_change = (model3$bic-model2$bic) / abs(model2$bic)
    
    #Only use model with 3 components if it substantially improves BIC
    if(fraction_change < 0.04){
      model = model2
    }else{
      model = model3
    }
    
    #Extract parameters of the two largest fitted components
    mean = model$parameters$mean
    mean = mean[(length(mean)-1):length(mean)]
    sigmasq = model$parameters$variance$sigmasq
    sigmasq = sigmasq[(length(sigmasq)-1):length(sigmasq)]
    params = data.frame(sample = sample_id, mean1 = mean[1], mean2 = mean[2], sigmasq1 = sigmasq[1], sigmasq2 = sigmasq[2])
    
    #Identify fraction of positive cells
    treshold = params$mean2-3*sqrt(params$sigmasq2)
    purity = length(which(stained[,channel] > treshold)) / nrow(stained)
    params = dplyr::mutate(params, channel = channel, purity = purity)
    
    results = rbind(results, params)
  }
  return(results)
}

changeFlowFrameNames <- function(flowframe, name_map){
  old_names = flowCore::colnames(flowframe)
  map_names = dplyr::left_join(data.frame(old_names = old_names, stringsAsFactors = FALSE), name_map, by = "old_names") %>% 
    dplyr::mutate(new_names = ifelse(is.na(new_names), old_names, new_names))
  flowCore::colnames(flowframe) = map_names$new_names
  return(flowframe)
}

addMetadataToFlowset <- function(flowset, metadata){
  
  #Add additional columns from metadata
  new_pData = dplyr::left_join(pData(flowset), metadata, by = "name")
  rownames(new_pData) = new_pData$name
  new_phenoData = as(new_pData, "AnnotatedDataFrame")
  
  #Update the varMetadata slot
  var_meta = varMetadata(new_phenoData)
  var_meta$labelDescription = rownames(var_meta)
  varMetadata(new_phenoData) = var_meta
  
  phenoData(flowset) = new_phenoData
  return(flowset)
}
