# Authors: Richard Schuster & Dan Wismer
#
# Date: September 1st, 2022
#
# Description: This script extracts national data to an area of interest (aoi)
#
# Inputs:  1. an area of interest (aoi) raster
#          2. an output folder location
#          3. required R libraries
#
# Outputs: 1. a 1km x 1km raster layer for each variable that intersects 
#             with the aoi
#          2. a csv that lists the species that intersect the aoi
#

# 1.0 Load packages ------------------------------------------------------------

## Start timer
start_time <- Sys.time()

## Package names
packages <- c("sf", "raster", "dplyr", "prioritizr", "devtools", "stringr")

## Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

## gdalUtils is now an archived package
if(!require(gdalUtils)){
  devtools:::install_github("gearslaboratory/gdalUtils")
}

## Load packages
library(sf)
library(raster)
library(dplyr)
library(gdalUtils) 
library(prioritizr)
library(stringr)
source("R/fct_matrix_intersect.R")
source("R/fct_matrix_to_raster.R")

# 2.0 Set up -------------------------------------------------------------------

## Set output folder ----
output_folder <- "data/output/Tiffs" # <--- CHANGE HERE FOR NEW PROJECT

## Set output csv ----
output_csv <- "data/csv/species.csv" # <--- CHANGE HERE FOR NEW PROJECT

## Read-in area of interest .tiff (aoi) ----
aoi_path <- "data/aoi/R1km_AOI.tif" # <--- CHANGE HERE FOR NEW AOI
aoi_1km <- raster(aoi_path) 
aoi_name <- names(aoi_1km)
aoi_1km0 <- aoi_1km 
aoi_1km0[aoi_1km0 == 1] <- 0

## Read-in national 1km grid (all of Canada, land only) ----
ncc_1km <- raster("data/national/_nccgrid/boundary.tif")
ncc_1km_idx <- ncc_1km
ncc_1km_idx[] <- 1:ncell(ncc_1km_idx) # 267,790,000 planning units

## Align aoi to same extent and same number of rows/cols as national grid ----
unlink("data/aoi/align/*", recursive = T, force = T) # clear folder
gdalUtils::align_rasters(unaligned = aoi_path,
                         reference = "data/national/_nccgrid/boundary.tif",
                         dstfile = paste0("data/aoi/align/", aoi_name, ".tif"))

## Get aligned AOI planning units ---- 
aoi_pu <- raster(paste0("data/aoi/align/", aoi_name, ".tif"))
# Create aoi_rij matrix: 9,880,260 planing units
aoi_rij <- prioritizr::rij_matrix(ncc_1km, stack(aoi_pu, ncc_1km_idx)) 
rownames(aoi_rij) <- c("AOI", "Idx")


# 3.0 Extract species to aoi ---------------------------------------------------

## ECCC Species at risk (theme) ----
natdata_rij <- readRDS("data/national/species/rij_SAR.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0, 
                 output_folder, "T_ECCC_SAR_", "INT1U")

## Nature Serve Canada Species at risk (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_NSC_SAR.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_NSC_SAR_", "INT1U")

## Nature Serve Canada Endemics (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_NSC_END.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_NSC_END_", "INT1U")

## Nature Serve Canada Common Species (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_NSC_SPP.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_NSC_SPP_", "INT1U")

## Amphibians (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_amph.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_IUCN_Amphibians_", "INT1U")

## Birds (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_bird.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_IUCN_Birds_", "INT1U")

## Eco regions (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_ecor.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_IUCN_EcoRegions_", "INT1U")

## Mammals (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_mamm.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_IUCN_Mammals_", "INT1U")

## Reptiles (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_rept.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_IUCN_Reptiles_", "INT1U")

# 4.0 Extract conservation variables to aoi ------------------------------------

## Carbon storage (weight) ----
natdata_r <- raster("data/national/carbon/Carbon_Mitchell_2021_t.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Carbon_storage")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "W_", "FLT4S")

## Carbon potential (weight) ----
natdata_r <- raster("data/national/carbon/Carbon_Potential_NFI_2011_CO2e_t_year.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Carbon_potential")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "W_", "FLT4S")

## Climate forward velocity (weight) ----
natdata_r <- raster("data/national/climate/fwdshortestpath.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Climate_forward_velocity")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "W_", "FLT4S")

## Climate refugia (weight) ----
natdata_r <- raster("data/national/climate/NA_combo_refugia_sum45.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Climate_refugia")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "W_", "FLT4S")

## Forest (theme) ----
natdata_r <- raster("data/national/forest/CA_forest_VLCE_2015_forest_only_ha_proj_scale.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Forest")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_LC_", "INT2U")

## Grassland (theme) ----
natdata_r <- raster("data/national/grassland/AAFC_LU2015_comb_masked_by_Prairie_grassland_comb.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Grassland")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_LC_", "INT2U")

## Human footprint (weight) ----
natdata_r <- raster("data/national/disturbance/cum_threat_int_proj2.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Human_disturbance")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "W_", "FLT4S")

## KBAs (weight) ----
natdata_r <- raster("data/national/kba/KBA.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Key_biodiversity_areas")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "W_", "FLT4S")

## Protected (include) ----
### canadian protected and conserved areas database (CPCAD)
CPCAD <- raster("data/national/protected/CPCAD.tif")
CPCAD[is.na(CPCAD[])] <- 0 # convert na values to 0
### NCC direct properties
NCC_direct <- raster("data/national/protected/NCC_direct.tif")
NCC_direct[is.na(NCC_direct[])] <- 0 
### NCC indirect properties
NCC_indirect <- raster("data/national/protected/NCC_indirect.tif") 
NCC_indirect[is.na(NCC_indirect[])] <- 0 
### raster math
CPCAD_NCC <- CPCAD + NCC_direct + NCC_indirect
CPCAD_NCC[CPCAD[] > 0] <- 1 # values > 0, convert to a value of 1

natdata_rij <- prioritizr::rij_matrix(ncc_1km, CPCAD_NCC)
rownames(natdata_rij) <- c("Protected")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "I_", "INT1U")
### remove objects to save RAM
rm(CPCAD, NCC_direct, NCC_indirect)
gc()

## Recreation (weight) ----
natdata_r <- raster("data/national/recreation/rec_pro_1a_norm.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Recreation")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "W_", "FLT4S")

## Freshwater (weight) ----
natdata_r <- raster("data/national/water/water_provision_2a_norm.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Freshwater")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "W_", "FLT4S")

## River length (theme) ----
natdata_r <- raster("data/national/water/grid_1km_water_linear_flow_length_1km.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("River_length")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_LC_", "FLT4S")

## Lakes (theme) ----
natdata_r <- raster("data/national/water/Lakes_CanVec_50k_ha.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Lakes")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_LC_", "FLT4S")

## Shoreline (theme) ----
natdata_r <- raster("data/national/water/Shoreline.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Shoreline_length")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_LC_", "FLT4S")

## Wetlands (theme) ----
natdata_r <- raster("data/national/wetlands/Wetland_comb_proj_diss_90m_Arc.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Wetlands")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
                 output_folder, "T_LC_", "FLT4S")

# Clear some RAM if possible
gc()

# 5.0 Generate species list csv ------------------------------------------------

## List files in output folder
output_tiffs <- list.files(output_folder, pattern='.tif$', full.names = T, 
                      recursive = F)

## Create empty species vectors to populate
ECCC_SAR <- c()
NSC_SAR <- c()
NSC_END <- c()
NSC_SPP <- c()
IUCN_Amphibians <- c()
IUCN_Birds <- c()
IUCN_Mammals <- c()
IUCN_Reptiles <- c()

## Populate species vectors ----
for (tiff in output_tiffs) {
  ### get file name
  name <- tools::file_path_sans_ext(basename(tiff))
  
  ### ECCC_SAR
  if(str_detect(name, "T_ECCC_SAR_")) {
    species_name <- unlist(str_split(name, "T_ECCC_SAR_"))[2]
    ECCC_SAR <- c(ECCC_SAR, species_name)
  }
  
  ### NSC_SAR
  if(str_detect(name, "T_NSC_SAR_")) {
    species_name <- unlist(str_split(name, "T_NSC_SAR_"))[2]
    NSC_SAR <- c(NSC_SAR, species_name)
  } 
  
  ### NSC_END
  if(str_detect(name, "T_NSC_END_")) {
    species_name <- unlist(str_split(name, "T_NSC_END_"))[2]
    NSC_END <- c(NSC_END, species_name)
  } 
  
  ### NSC_SPP
  if(str_detect(name, "T_NSC_SPP_")) {
    species_name <- unlist(str_split(name, "T_NSC_SPP_"))[2]
    NSC_SPP <- c(NSC_SPP, species_name)
  } 
  
  ### IUCN_Amphibians
  if(str_detect(name, "T_IUCN_Amphibians_")) {
    species_name <- unlist(str_split(name, "T_IUCN_Amphibians_"))[2]
    IUCN_Amphibians <- c(IUCN_Amphibians, species_name)
  } 
  
  ### IUCN_Birds
  if(str_detect(name, "T_IUCN_Birds_")) {
    species_name <- unlist(str_split(name, "T_IUCN_Birds_"))[2]
    IUCN_Birds <- c(IUCN_Birds, species_name)
  }
  
  ### IUCN_Mammals
  if(str_detect(name, "T_IUCN_Mammals_")) {
    species_name <- unlist(str_split(name, "T_IUCN_Mammals_"))[2]
    IUCN_Mammals <- c(IUCN_Mammals, species_name)
  }  
  
  ### IUCN_Reptiles
  if(str_detect(name, "T_IUCN_Reptiles_")) {
    species_name <- unlist(str_split(name, "T_IUCN_Reptiles_"))[2]
    IUCN_Reptiles <- c(IUCN_Reptiles, species_name)
  }  
}


## Create csv ----
## Get the length of the longest vector 
max_ln <- max(c(length(ECCC_SAR), length(NSC_END),
              length(NSC_END), length(NSC_SPP),
              length(IUCN_Amphibians), length(IUCN_Birds),
              length(IUCN_Mammals), length(IUCN_Reptiles)))

## Build data.frame 
species <- data.frame(ECC_SAR = c(ECCC_SAR, rep(NA, max_ln - length(ECCC_SAR))),
                      NSC_END = c(NSC_END, rep(NA, max_ln - length(NSC_END))),
                      NSC_SAR = c(NSC_SAR, rep(NA, max_ln - length(NSC_SAR))),
                      NSC_SPP = c(NSC_SPP, rep(NA, max_ln - length(NSC_SPP))),
                      IUCN_Amphibians = c(IUCN_Amphibians, rep(NA, max_ln - length(IUCN_Amphibians))),
                      IUCN_Birds = c(IUCN_Birds, rep(NA, max_ln - length(IUCN_Birds))),
                      IUCN_Mammals = c(IUCN_Mammals, rep(NA, max_ln - length(IUCN_Mammals))),
                      IUCN_Reptiles = c(IUCN_Reptiles, rep(NA, max_ln - length(IUCN_Reptiles))))

## Write species data.frame to disk 
write.csv(species, 
          output_csv, 
          row.names = FALSE)

# 6.0 Clear R environment ------------------------------------------------------ 

## End timer
end_time <- Sys.time()
end_time - start_time

rm(list=ls())
gc()
