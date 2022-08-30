#' Writes national data rij matrix to raster.
#' 
#' @description
#' The matrix_to_raster function takes the national-AOI rij matrix  and 
#' populates a template raster by setting values based on the cell index. 
#' 
#' @param ncc_1km_idx a [raster] that has the NCC national grid index as values.
#' 
#' @param natdata_intersect a sparse matrix of class [dgCMatrix]. The 
#' natdata-aoi intersect matrix.
#' 
#' @param aoi_1km0 an aoi [raster]. All values must be set to 0. This raster is 
#' used to crop the national-aoi raster back to the original aoi and also create
#' the binary output where aoi values == 0 and national data == 1.
#' 
#' @param output_folder a [character] path to output raster layers.
#' 
#' @param prefix a [character] string to append before the raster name. Helpful
#' for keeping standard naming. Ex. T_ECC_SAR_. Where T == theme, 
#' ECC_SAR == source. 

matrix_to_raster = function(ncc_1km_idx, natdata_intersect, aoi_1km0, 
                            output_folder, prefix) {
  
  # Set up placeholder raster
  natdata_raster <- ncc_1km_idx
  natdata_rasters <- list()
  
  # Set up counter for print message
  len =  (nrow(natdata_intersect)-2)
  counter = 1
  
  # Loop through matrix, exclude AOI and Idx rows
  for (i in 1:(nrow(natdata_intersect)-2)) {
    natdata_raster[] <- NA
    name <- rownames(natdata_intersect)[i]
    print(paste0("... ", counter, " of ", len, ": ",  name))
    natdata_raster[natdata_intersect["Idx",]] <- natdata_intersect[i,]
    names(natdata_raster) <- name
    ## crop raster to AOI
    cropped <- raster::crop(natdata_raster, aoi_1km0)
    ## create binary raster (aoi == 0, natdata = 1) 
    feature <- raster::mosaic(cropped, aoi_1km0, fun = "max")
    writeRaster(feature, paste0(output_folder, "/", prefix, name,".tif"), overwrite = TRUE, datatype = "INT2S")
    counter = counter + 1
  }
}
