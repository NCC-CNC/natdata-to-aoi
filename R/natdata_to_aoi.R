# Authors: Richard Schuster & Dan Wismer
# Date: August 26th, 2022
# Description: extract national data to an area of interest (aoi) 

library(sf)
library(raster)
library(dplyr)
library(gdalUtils)
library(prioritizr)
source("R/fct_matrix_intersect.R")
source("R/fct_matrix_to_raster.R")

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

# Align aoi to same extent and same number of rows/cols as national grid ----
gdalUtils::align_rasters(unaligned = aoi_path,
                         reference = "data/national/boundary.tif",
                         dstfile = paste0("data/aoi/align/", aoi_name, ".tif"))

# Get aligned AOI planning units ---- 
aoi_pu <- raster(paste0("data/aoi/align/", aoi_name, ".tif"))
# Create rij matrix
aoi_rij <- prioritizr::rij_matrix(ncc_1km, stack(aoi_pu, ncc_1km_idx)) 
rownames(aoi_rij) <- c("AOI", "Idx")

# Convert sparse matrix back to raster ----

## ECCC Species at risk ----
ECC_SAR <- readRDS("data/national/rij_SAR.rds")
ECC_SAR_intersect <- matrix_intersect(ECC_SAR, aoi_rij) 
matrix_to_raster(ncc_1km_idx, ECC_SAR_intersect, output_folder, "T_ECC_SAR_")
## clean up to save RAM
rm(ECC_SAR_intersect)

## Nature Serve Canada Species at risk ----
NSC_SAR <- readRDS( "data/national/rij_NSC_SAR.rds")
NSC_SAR_intersect <- matrix_intersect(NSC_SAR, aoi_rij) 
matrix_to_raster(ncc_1km_idx, NSC_SAR_intersect, output_folder, "T_NSC_SAR_")
## clean up to save RAM
rm(NSC_SAR_intersect)

## Nature Serve Canada Endemics ----
NSC_END <- readRDS( "data/national/rij_NSC_END.rds")
NSC_END_intersect <- matrix_intersect(NSC_END, aoi_rij) 
matrix_to_raster(ncc_1km_idx, NSC_END_intersect, output_folder, "T_NSC_END_")
## clean up to save RAM
rm(NSC_END_intersect)

## Nature Serve Canada Common Species ----
NSC_SPP <- readRDS( "data/national/rij_NSC_SPP.rds")
NSC_SPP_intersect <- matrix_intersect(NSC_SPP, aoi_rij) 
matrix_to_raster(ncc_1km_idx, NSC_SPP_intersect, output_folder, "T_NSC_SPP_")
## clean up to save RAM
rm(NSC_SPP_intersect)

#-------------------------------------------------------------------------------

# Generate species list csv ----

## NCC planning units
ncc_pu <- ncc_1km[]
## AOI planning units ----
aoi_pu <- raster(paste0("data/aoi/align/", aoi_name, ".tif")) %>%
  getValues()  %>% .[!is.na(ncc_pu)]
aoi_pu <- ifelse(!is.na(aoi_pu), 1, 0)

## ECC Species at Risk matrix multiplication ----
ECC_SAR_X_aoi_pu <- ECC_SAR %*% aoi_pu
ECC_SAR_aoi_pu <- names(ECC_SAR_X_aoi_pu[rowSums(ECC_SAR_X_aoi_pu) > 0, ])

## Nature Serve Canada Endemics matrix multiplication ----
NSC_END_X_aoi_pu <- NSC_END %*% aoi_pu
NSC_END_aoi_pu <- names(NSC_END_X_aoi_pu[rowSums(NSC_END_X_aoi_pu) > 0, ])

## Nature Serve Canada Species at Risk matrix multiplication ----
NSC_SAR_X_aoi_pu <- NSC_SAR %*% aoi_pu
NSC_SAR_aoi_pu <- names(NSC_SAR_X_aoi_pu[rowSums(NSC_SAR_X_aoi_pu) > 0, ])

## Nature Serve Canada Common Species matrix multiplication ----
NSC_SPP_X_aoi_pu <- NSC_SPP %*% aoi_pu
NSC_SPP_aoi_pu <- names(NSC_SPP_X_aoi_pu[rowSums(NSC_SPP_X_aoi_pu) > 0, ])

## Create csv ----
## Get the length of the longest vector 
max_ln <- max(c(length(ECC_SAR_aoi_pu ), length(NSC_END_aoi_pu )),
              (length(NSC_END_aoi_pu)), length(NSC_SPP_aoi_pu))

## Build data.frame 
species <- data.frame(ECC_SAR = c(ECC_SAR_aoi_pu,rep(NA, max_ln - length(ECC_SAR_aoi_pu))),
                      NSC_END = c(NSC_END_aoi_pu,rep(NA, max_ln - length(NSC_END_aoi_pu))),
                      NSC_SAR = c(NSC_SAR_aoi_pu,rep(NA, max_ln - length(NSC_SAR_aoi_pu))),
                      NSC_SPP = c(NSC_SPP_aoi_pu,rep(NA, max_ln - length(NSC_SPP_aoi_pu))))

## Write species data.frame to disk ----
write.csv(species, 
          paste0("data/csv/", aoi_name, "_species.csv"), 
          row.names = FALSE) 

