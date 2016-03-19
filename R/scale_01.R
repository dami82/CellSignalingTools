scale_01 <-
function(num_matrix) {
  my_min <- min(num_matrix)
  my_max <- max(num_matrix)
  res_mat <- (num_matrix - my_min) / (my_max - my_min)
  return(res_mat)
}
