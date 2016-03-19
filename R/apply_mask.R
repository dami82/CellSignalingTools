apply_mask <-
function(green_image, dapi_mask, bckgrnd = 0.0){
  if(ncol(green_image) == ncol(dapi_mask) & 
     nrow(green_image) == nrow(dapi_mask) &
     is.numeric(bckgrnd) &
     bckgrnd < 0.9 &
     bckgrnd >= 0){
    #calculate_dim
    dapi_dims <- get_nucl_dims(dapi_range = dapi_mask)
    #remove_dust
    my_mask <- remove_dust(dapi_mask,nuclei_dims = dapi_dims)  
    #
    message("applying mask...")
    #split green image
    green_img <- scale_01(green_image)
    green_img[green_img < bckgrnd] <- 0 
    #
    nucl_signal <- sapply(1:ncol(my_mask), (function(xx){
      sapply(1:nrow(my_mask), (function(yy){
        if(my_mask[yy,xx] == TRUE) {
          green_img[yy,xx]    
        } else {
          FALSE
        }
      }))
    }))
    cyto_signal <- sapply(1:ncol(my_mask), (function(xx){
      sapply(1:nrow(my_mask), (function(yy){
        if(my_mask[yy,xx] == FALSE) {
          green_img[yy,xx]    
        } else {
          FALSE
        }
      }))
    }))
    #
    result <- list()
    result[["nucl"]] <- nucl_signal
    result[["cyto"]] <- cyto_signal
    result[["input"]] <- green_img 
    return(result)
  } else {
    message("Bad input!")
  }  
}
