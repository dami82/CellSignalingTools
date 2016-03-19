get_nucl_dims <-
function(dapi_range) {
  #init
  #row by row
  max_width <- max(sapply(1:nrow(dapi_range), (function(i){
    max_x <- list()
    tmp_vct <- which(dapi_range[i,] == TRUE)
    if (length(tmp_vct)>1) {
      k = 1
      for(j in 2:length(tmp_vct)){
        if(tmp_vct[j] - tmp_vct[j-1] == 1) {
          max_x[[k]] <- 1
        } else {
          k <- 1
          max_x[[k]] <- 1
        }
        k = k+1  
      }
    }
    length(max_x)
  })))
  #col by col
  max_height <- max(sapply(1:ncol(dapi_range), (function(i){
    max_y <- list()
    tmp_vct <- which(dapi_range[,i] == TRUE)
    if (length(tmp_vct)>1) {
      k = 1
      for(j in 2:length(tmp_vct)){
        if(tmp_vct[j] - tmp_vct[j-1] == 1) {
          max_y[[k]] <- 1
        } else {
          k <- 1
          max_y[[k]] <- 1
        }
        k = k+1  
      }
    }
    length(max_y)
  })))
  return(c(max_width, max_height))
}
