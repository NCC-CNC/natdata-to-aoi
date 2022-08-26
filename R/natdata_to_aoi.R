#
# Authors: Richard Schuster & Dan Wismer
# Date: August 26th, 2022
# Description:
#
#
#
#

# Import libraries and user functions
library(sf)
library(raster)
source("R/fct_res_processing.R")

# Read-in national 1km grid (all of Canada) ----
ncc_1km <- raster("data/national/Constant_1km.tif")
ncc_pu <- ncc_1km[]
ncc_pu <- ncc_pu[!is.na(ncc_pu)]


# Read-in area of interest .tiff (aoi) ----
aoi_pu <- raster("data/aoi/R1km_AOI.tif") %>%
  getValues() %>%
  .[!is.na(ncc_pu)] 
# Reclass na to 0
aoi_pu <- ifelse(!is.na(aoi_pu), 1, 0)

# Read-in national data ----
SAR <- readRDS("data/national/rij_SAR.rds") 
# amph <- readRDS( "data/national/rij_amph.rds") 
# bird <- readRDS( "data/national/rij_bird.rds") 
# mamm <- readRDS( "data/national/rij_mamm.rds") 
# rept <- readRDS( "data/national/rij_rept.rds") 
# ecor <- readRDS( "data/national/rij_ecor.rds") 
# NSC_END <- readRDS( "data/national/rij_NSC_END.rds") 
# NSC_SAR <- readRDS( "data/national/rij_NSC_SAR.rds") 
# NSC_SPP <- readRDS( "data/national/rij_NSC_SPP.rds")

# Matrix multiplication 
pu_SAR <- SAR %*% aoi_pu 
