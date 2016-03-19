remove_dust <-
function(signal_map, nuclei_dims) {
  d_param <- 10
  min_val <- 5
  #
  gap_x <- as.integer(nuclei_dims[1]/d_param)
  gap_y <- as.integer(nuclei_dims[2]/d_param)
  #
  my_map <- signal_map
  #
  message("adjusting mask... this may take a couple of minutes...")
  nu_map <- sapply(1:(ncol(my_map)), (function(cl){
    sapply(1:(nrow(my_map)), (function(rw){
      #
      top_range <- rw - gap_y
      if (top_range < 1) {top_range <- 1 }
      bot_range <- rw + gap_y
      if (bot_range > nrow(my_map)) {bot_range <- nrow(my_map) }
      lef_range <- cl - gap_x
      if(lef_range < 1) {lef_range <- 1 }  
      rig_range <- cl + gap_x
      if(rig_range > ncol(my_map)) {rig_range <- ncol(my_map) }
      #
      if(my_map[rw,cl] == FALSE) {
        #
        if ((sum(my_map[rw,lef_range:cl]) > min_val &&
             sum(my_map[rw,cl:rig_range]) > min_val) | 
            (sum(my_map[top_range:rw,cl]) > min_val &&
             sum(my_map[rw:bot_range,cl]) > min_val)){
          TRUE
        } else {
          FALSE
        }
      } else {
        TRUE
      }
    }))
  }))
  #
  my_map <- nu_map
  nu_map <- sapply(1:(ncol(my_map)), (function(cl){
    sapply(1:(nrow(my_map)), (function(rw){
      #
      top_range <- rw - gap_y
      if (top_range < 1) {top_range <- 1 }
      bot_range <- rw + gap_y
      if (bot_range > nrow(my_map)) {bot_range <- nrow(my_map) }
      lef_range <- cl - gap_x
      if(lef_range < 1) {lef_range <- 1 }  
      rig_range <- cl + gap_x
      if(rig_range > ncol(my_map)) {rig_range <- ncol(my_map) }
      #
      if(my_map[rw,cl] == TRUE) {
        #
        if ((sum(my_map[rw,lef_range:rig_range]) < gap_x) | 
            (sum(my_map[top_range:bot_range,cl]) < gap_y)){
          FALSE
        } else {
          TRUE
        }
      } else {
        FALSE
      }
    }))
  }))
  return(nu_map)
}
