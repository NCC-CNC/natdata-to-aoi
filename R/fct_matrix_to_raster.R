matrix_to_raster = function(ncc_1km_idx, natdata_intersect, output_folder, prefix) {
  
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