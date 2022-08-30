#' Finds national data records that overlap AOI
#' 
#' @description
#' The matrix_intersect function uses indexing to subset rows that are greater
#' than 0.
#' 
#' @param natdata_rij a sparse matrix of class [dgCMatrix]. National data,
#' example `rij_NSC_SAR.rds`
#' 
#' @param aoi_rij a sparse matrix of class [dgCMatrix]. Area of interest 
#' `rij_matrix`

matrix_intersect <- function(natdata_rij, aoi_rij) {
  # append AOI and Idx to SAR matrix
  natdata_aoi_ncc <- rbind(natdata_rij, aoi_rij) 
  # get all values greater than (gt) 0 
  natdata_aoi_ncc_gt0 <- natdata_aoi_ncc[, natdata_aoi_ncc["AOI",] > 0]
  # subset species that `intersect` aoi 
  natdata_intersect <- natdata_aoi_ncc_gt0[rowSums(natdata_aoi_ncc_gt0) > 0, ] 
}