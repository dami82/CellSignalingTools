get_auc <-
function(signal_matrix){
  signal_mat <- signal_matrix
  colnames(signal_mat) <- c("pos","signal")
  #
  auc <- sum(sapply(2:nrow(signal_mat), (function(x){
    mean_y <- mean(signal_mat[x,"signal"],signal_mat[x-1,"signal"])
    x1 <- signal_mat[x,"pos"]
    x2 <- signal_mat[x-1,"pos"]
    x_base <- max(x1, x2) - min(x1, x2)
    auc <- mean_y * x_base       
  })))
  #
  return(auc)
}
