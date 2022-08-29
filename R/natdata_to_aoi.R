# Authors: Richard Schuster & Dan Wismer
# Date: August 26th, 2022
# Description: extacts national data to an area of intest (aoi) 

library(sf)
library(raster)
library(dplyr)
library(gdalUtils)
library(prioritizr)

# Set output folder
output_folder <- "data/output/Tiffs" # <---------- CHANGE HERE FOR NEW PROJECT

# Read-in area of interest .tiff (aoi) ----
aoi_path <- "data/aoi/R1km_AOI.tif" # <---------------- CHANGE HERE FOR NEW AOI
aoi_1km <- raster(aoi_path) 
aoi_name <- names(aoi_1km)
aoi_1km0 <- aoi_1km 
aoi_1km0[aoi_1km0 == 1] <- 0

# Read-in national 1km grid (all of Canada) ----
ncc_1km <- raster("data/national/boundary.tif")
ncc_1km_idx <- ncc_1km
ncc_1km_idx[] <- 1:ncell(ncc_1km_idx) # could restrict to non-NA's is desired

# NCC planning units
ncc_pu <- ncc_1km[]

# Align aoi to same extent and same number of rows/cols as national grid ----
gdalUtils::align_rasters(unaligned = aoi_path,
                         reference = "data/national/boundary.tif",
                         dstfile = paste0("data/aoi/align/", aoi_name, ".tif"))

# Get aligned AOI planning units ---- 
aoi_pu <- raster(paste0("data/aoi/align/", aoi_name, ".tif"))
# Create rij matrix
aoi_rij <- prioritizr::rij_matrix(ncc_1km, stack(aoi_pu, ncc_1km_idx)) 
rownames(aoi_rij) <- c("AOI", "Idx")

# Read-in national data ----
ECC_SAR <- readRDS("data/national/rij_SAR.rds")
NSC_END <- readRDS( "data/national/rij_NSC_END.rds")
NSC_SAR <- readRDS( "data/national/rij_NSC_SAR.rds")
NSC_SPP <- readRDS( "data/national/rij_NSC_SPP.rds")

## ECCC Species at risk ----
SAR_aoi_ncc <- rbind(ECC_SAR, aoi_rij) # append AOI and Idx to SAR matrix
SAR_aoi_ncc_gt0 <- SAR_aoi_ncc[, SAR_aoi_ncc["AOI",] > 0] # get all values greater than (gt) 0 (ECCC SAR output: 486 x 12990)
SAR_intersect <- SAR_aoi_ncc_gt0[rowSums(SAR_aoi_ncc_gt0) > 0,] # subset species that `intersect` aoi 

# Convert sparse matrix back to raster ----
natdata_raster <- ncc_1km_idx
natdata_rasters <- list()

## Loop through matrix, exclude AOI and Idx rows
len =  (nrow(SAR_intersect)-2)
counter = 1
for (i in 1:(nrow(SAR_intersect)-2)) {
  natdata_raster[] <- NA
  name <- rownames(SAR_intersect)[i]
  print(paste0("... ", counter, " of ", len, ": ",  name))
  natdata_raster[SAR_intersect["Idx",]] <- SAR_intersect[i,]
  names(natdata_raster) <- name
  cropped <- raster::crop(natdata_raster, aoi_1km0)
  feature <- raster::mosaic(cropped, aoi_1km0, fun = "max")
  writeRaster(feature, paste0(output_folder, "/", name,".tif"), overwrite = TRUE, datatype = "INT2S")
  counter = counter + 1
}

#-------------------------------------------------------------------------------

# Build out species list ----
aoi_pu <- raster(paste0("data/aoi/align/", aoi_name, ".tif")) %>%
  getValues()  %>% .[!is.na(ncc_pu)]
# Reclass na to 0
aoi_pu <- ifelse(!is.na(aoi_pu), 1, 0)

## ECC Species at Risk
ECC_SAR_X_aoi_pu <- ECC_SAR %*% aoi_pu
ECC_SAR_aoi_pu <- names(ECC_SAR_X_aoi_pu[rowSums(ECC_SAR_X_aoi_pu) > 0, ])

## Nature Serve Canada Endemics
NSC_END_X_aoi_pu <- NSC_END %*% aoi_pu
NSC_END_aoi_pu <- names(NSC_END_X_aoi_pu[rowSums(NSC_END_X_aoi_pu) > 0, ])

## Nature Serve Canada Species at Risk
NSC_SAR_X_aoi_pu <- NSC_SAR %*% aoi_pu
NSC_SAR_aoi_pu <- names(NSC_SAR_X_aoi_pu[rowSums(NSC_SAR_X_aoi_pu) > 0, ])

## Nature Serve Canada Common Species
NSC_SPP_X_aoi_pu <- NSC_SPP %*% aoi_pu
NSC_SPP_aoi_pu <- names(NSC_SPP_X_aoi_pu[rowSums(NSC_SPP_X_aoi_pu) > 0, ])

# Create csv ----
## Get the length of the longest vector ----
max_ln <- max(c(length(ECC_SAR_aoi_pu ), length(NSC_END_aoi_pu )),
              (length(NSC_END_aoi_pu)), length(NSC_SPP_aoi_pu))

## Build data.frame ----
species <- data.frame(ECC_SAR = c(ECC_SAR_aoi_pu,rep(NA, max_ln - length(ECC_SAR_aoi_pu))),
                      NSC_END = c(NSC_END_aoi_pu,rep(NA, max_ln - length(NSC_END_aoi_pu))),
                      NSC_SAR = c(NSC_SAR_aoi_pu,rep(NA, max_ln - length(NSC_SAR_aoi_pu))),
                      NSC_SPP = c(NSC_SPP_aoi_pu,rep(NA, max_ln - length(NSC_SPP_aoi_pu))))

## Write species data.frame to disk ----
write.csv(species, 
          paste0("data/csv/", aoi_name, "_species.csv"), 
          row.names = FALSE) 

