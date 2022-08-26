#
# Authors: Richard Schuster & Dan Wismer
# Date: August 26th, 2022
# Description: remember to fill this in
#
#

# Import libraries and user functions
library(sf)
library(raster)
library(dplyr)
library(gdalUtils)
source("R/fct_res_processing.R")

# Read-in national 1km grid (all of Canada) ----
ncc_1km <- raster("data/national/boundary.tif")
ncc_pu <- ncc_1km[]

# Read-in area of interest .tiff (aoi) ----
aoi_1km <- raster("data/aoi/R1km_AOI.tif")

# Align aoi to same extent and same number of rows/cols as national grid ----
gdalUtils::align_rasters(unaligned = "data/aoi/R1km_AOI.tif",
                         reference = "data/national/boundary.tif",
                         dstfile = "data/aoi/align/R1km_AOI_Aligned.tif")

# Get numeric vector of aoi planning units 
aoi_pu <- raster("data/aoi/align/R1km_AOI_Aligned.tif") %>%
  getValues()  %>% .[!is.na(ncc_pu)] 
# Reclass na to 0
aoi_pu <- ifelse(!is.na(aoi_pu), 1, 0)

# Read-in national data ----
SAR <- readRDS("data/national/rij_SAR.rds") 
NSC_END <- readRDS( "data/national/rij_NSC_END.rds")
NSC_SAR <- readRDS( "data/national/rij_NSC_SAR.rds")
NSC_SPP <- readRDS( "data/national/rij_NSC_SPP.rds")

# Matrix multiplication 
SAR_X_aoi_pu <- SAR %*% aoi_pu
SAR_aoi_pu <- SAR_X_aoi_pu[rowSums(SAR_X_aoi_pu) > 0, ]
