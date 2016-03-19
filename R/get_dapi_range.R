get_dapi_range <-
function(dapi_measures, cutoff = "auto"){
  #start
  dapi <- scale_01(dapi_measures)
  init_cutoff <- 0.4
  med_low <- mean(dapi[dapi<=init_cutoff])
  sd_low <- sd(dapi[dapi<=init_cutoff])
  final_cutoff <- med_low + 2 * sd_low
  #
  if (is.numeric(cutoff)){
    if (cutoff > 0 && cutoff < 1) {
      final_cutoff <- cutoff  
    }
  }
  #
  if(is.vector(dapi_measures)){
    dapi_init <- min(which(dapi>final_cutoff))
    dapi_stop <- max(which(dapi>final_cutoff))
    return (c(1,dapi_init, dapi_stop, length(dapi)))  
  } else {
    dapi_plus <- dapi>final_cutoff
    #remove_dust
    return(dapi_plus)  
  }
}
