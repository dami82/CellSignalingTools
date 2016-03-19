nc_ratio_tiff <-
function(signal_data, signal_cutoff = 0.0){
  if (is.list(signal_data)) {
    #
    plot_signal(signal_data)
    #
    nucl_vector <- as.vector(signal_data$nucl)
    tot_nucl <- sum(nucl_vector[nucl_vector>signal_cutoff])
    #
    cyto_vector <- as.vector(as.vector(signal_data$cyto))
    tot_cyto <- sum(cyto_vector[cyto_vector>signal_cutoff])
    #
    image_ratio <- tot_nucl / tot_cyto
    return(image_ratio)
  }
}
