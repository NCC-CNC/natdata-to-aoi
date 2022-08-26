res_processing <- function(rij_in = NULL, pu = NULL) {
  # matrix multiplication
  tmp <- rij_in %*% pu 
  # get column names
  nms <- tmp@Dimnames[[1]]
  # create named vector
  tmp <- as.vector(tmp)
  names(tmp) <- nms
  return(tmp)
}