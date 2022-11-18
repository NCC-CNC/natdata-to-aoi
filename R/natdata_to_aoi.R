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
packages <- c("sf", "raster", "dplyr", "prioritizr", 
              "sp", "stringr", "gdalUtilities")

## Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

## Load packages
library(sf)
library(raster)
library(dplyr)
library(prioritizr)
library(sp)
library(stringr)
library(gdalUtilities) 
source("R/fct_matrix_intersect.R")
source("R/fct_matrix_to_raster.R")

# 2.0 Set up -------------------------------------------------------------------

## Set output folder and AOI ----
root_folder <- "data/output" # <--- CHANGE HERE FOR NEW PROJECT
aoi_path <- "data/aoi/AOI.tif" # <--- CHANGE HERE FOR NEW AOI

## Create output folder directory ----
dir.create(file.path(root_folder, "Variables"))
dir.create(file.path(root_folder, "Variables", "Excludes"))
dir.create(file.path(root_folder, "Variables", "Includes"))
dir.create(file.path(root_folder, "Variables", "Themes"))
dir.create(file.path(root_folder, "Variables", "Weights"))

dir.create(file.path(root_folder, "Variables", "Themes", "ECCC_CH"))
dir.create(file.path(root_folder, "Variables", "Themes", "ECCC_SAR"))
dir.create(file.path(root_folder, "Variables", "Themes", "IUCN_AMPH"))
dir.create(file.path(root_folder, "Variables", "Themes", "IUCN_BIRD"))
dir.create(file.path(root_folder, "Variables", "Themes", "IUCN_MAMM"))
dir.create(file.path(root_folder, "Variables", "Themes", "IUCN_REPT"))
dir.create(file.path(root_folder, "Variables", "Themes", "LC"))
dir.create(file.path(root_folder, "Variables", "Themes", "NSC_END"))
dir.create(file.path(root_folder, "Variables", "Themes", "NSC_SAR"))
dir.create(file.path(root_folder, "Variables", "Themes", "NSC_SPP"))

## Read-in area of interest .tiff (aoi) ----
aoi_1km <- raster(aoi_path) 
aoi_name <- names(aoi_1km)
aoi_1km0 <- aoi_1km 
aoi_1km0[aoi_1km0 == 1] <- 0

## Read-in national 1km grid (all of Canada) ----
ncc_1km <- raster("data/national/_nccgrid/NCC_PU.tif")
ncc_1km_idx <- ncc_1km
ncc_1km_idx[] <- 1:ncell(ncc_1km_idx) # 267,790,000 available planning units

## Align aoi to same extent and same number of rows/cols as national grid ----
unlink("data/aoi/align/*", recursive = T, force = T) # clear folder
### get spatial properties of ncc grid
proj4_string <- sp::proj4string(ncc_1km) # projection string
bbox <- raster::bbox(ncc_1km) # bounding box
### variables for gdalwarp
te <- c(bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2]) # xmin, ymin, xmax, ymax
ts <- c(raster::ncol(ncc_1km), raster::nrow(ncc_1km)) # ncc grid: columns/rows
### gdalUtilities::gdalwarp does not require a local GDAL installation ----
gdalUtilities::gdalwarp(srcfile = aoi_path,
                        dstfile = paste0("data/aoi/align/", aoi_name, ".tif"),
                        te = te,
                        t_srs = proj4_string,
                        ts = ts,
                        overwrite = TRUE)

## Get aligned AOI planning units ---- 
aoi_pu <- raster(paste0("data/aoi/align/", aoi_name, ".tif"))
# Create aoi_rij matrix: 11,010,932 planing units activated 
aoi_rij <- prioritizr::rij_matrix(ncc_1km, stack(aoi_pu, ncc_1km_idx)) 
rownames(aoi_rij) <- c("AOI", "Idx")

# 3.0 Extract species to aoi ---------------------------------------------------

## ECCC Critical Habitat (theme) ----
natdata_rij <- readRDS("data/national/species/rij_ECCC_CH.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0, 
  paste0(root_folder, "/Variables/Themes/ECCC_CH"), "T_ECCC_CH_", "INT1U")

## ECCC Species at risk (theme) ----
natdata_rij <- readRDS("data/national/species/rij_ECCC_SAR.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0, 
  paste0(root_folder, "/Variables/Themes/ECCC_SAR"), "T_ECCC_SAR_", "INT1U")

## IUCN Amphibians (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_IUCN_AMPH.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/IUCN_AMPH"), "T_IUCN_AMPH_", "INT1U")

## IUCN Birds (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_IUCN_BIRD.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/IUCN_BIRD"), "T_IUCN_BIRD_", "INT1U")

## IUCN Mammals (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_IUCN_MAMM.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/IUCN_MAMM"), "T_IUCN_MAMM_", "INT1U")

## IUCN Reptiles (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_IUCN_REPT.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/IUCN_REPT"), "T_IUCN_REPT_", "INT1U")

## Nature Serve Canada Endemics (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_NSC_END.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/NSC_END"), "T_NSC_END_", "INT1U")

## Nature Serve Canada Species at risk (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_NSC_SAR.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/NSC_SAR"), "T_NSC_SAR_", "INT1U")

## Nature Serve Canada Common Species (theme) ----
natdata_rij <- readRDS( "data/national/species/rij_NSC_SPP.rds")
matrix_overlap <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/NSC_SPP"), "T_NSC_SPP_", "INT1U")

# 4.0 Extract conservation variables to aoi ------------------------------------

## Forest (theme) ----
natdata_r <- raster("data/national/forest/CA_forest_VLCE_2015_forest_only_ha_proj_scale.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Forest")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/LC"), "T_LC_", "INT2U")

## Grassland (theme) ----
natdata_r <- raster("data/national/grassland/AAFC_LU2015_comb_masked_by_Prairie_grassland_comb.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Grassland")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/LC"), "T_LC_", "INT2U")

## River length (theme) ----
natdata_r <- raster("data/national/water/grid_1km_water_linear_flow_length_1km.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("River_length")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/LC"), "T_LC_", "FLT4S")

## Lakes (theme) ----
natdata_r <- raster("data/national/water/Lakes_CanVec_50k_ha.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Lakes")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/LC"), "T_LC_", "FLT4S")

## Shoreline (theme) ----
natdata_r <- raster("data/national/water/Shoreline.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Shoreline_length")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/LC"), "T_LC_", "FLT4S")

## Wetlands (theme) ----
natdata_r <- raster("data/national/wetlands/Wetland_comb_proj_diss_90m_Arc.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Wetlands")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Themes/LC"), "T_LC_", "FLT4S")

## Carbon storage (weight) ----
natdata_r <- raster("data/national/carbon/Carbon_Mitchell_2021_t.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Carbon_storage")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Weights"), "W_", "FLT4S")

## Carbon potential (weight) ----
natdata_r <- raster("data/national/carbon/Carbon_Potential_NFI_2011_CO2e_t_year.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Carbon_potential")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Weights"), "W_", "FLT4S")

## Climate forward velocity (weight) ----
natdata_r <- raster("data/national/climate/fwdshortestpath.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Climate_forward_velocity")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Weights"), "W_", "FLT4S")

## Climate refugia (weight) ----
natdata_r <- raster("data/national/climate/NA_combo_refugia_sum45.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Climate_refugia")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Weights"), "W_", "FLT4S")

## Human footprint (weight) ----
natdata_r <- raster("data/national/disturbance/CDN_HF_cum_threat_20221031_NoData.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Human_footprint")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Weights"), "W_", "FLT4S")

## KBAs (weight) ----
natdata_r <- raster("data/national/kba/KBA.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Key_biodiversity_areas")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Weights"), "W_", "FLT4S")

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
  paste0(root_folder, "/Variables/Includes"), "I_", "INT1U")
### remove objects to save RAM
rm(CPCAD, NCC_direct, NCC_indirect)
gc()

## Recreation (weight) ----
natdata_r <- raster("data/national/recreation/rec_pro_1a_norm.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Recreation")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Weights"), "W_", "FLT4S")

## Freshwater (weight) ----
natdata_r <- raster("data/national/water/water_provision_2a_norm.tif")
natdata_rij <- prioritizr::rij_matrix(ncc_1km, natdata_r)
rownames(natdata_rij) <- c("Freshwater")
matrix_overlap  <- matrix_intersect(natdata_rij, aoi_rij) 
matrix_to_raster(ncc_1km_idx, matrix_overlap, aoi_1km0,
  paste0(root_folder, "/Variables/Weights"), "W_", "FLT4S")

# Clear some RAM if possible
gc()

# 5.0 Generate species list csv ------------------------------------------------

## List files in root folder
output_tiffs <- list.files(root_folder, pattern='.tif$', full.names = T, 
                      recursive = T)

## Create empty species vectors to populate
ECCC_SAR <- c()
ECCC_CH <- c()
NSC_SAR <- c()
NSC_END <- c()
NSC_SPP <- c()
IUCN_AMPH <- c()
IUCN_BIRD <- c()
IUCN_MAMM <- c()
IUCN_REPT <- c()

## Populate species vectors ----
for (tiff in output_tiffs) {
  ### get file name
  name <- tools::file_path_sans_ext(basename(tiff))
  
  ### ECCC_SAR
  if(str_detect(name, "T_ECCC_SAR_")) {
    species_name <- unlist(str_split(name, "T_ECCC_SAR_"))[2]
    ECCC_SAR <- c(ECCC_SAR, species_name)
  }
  
  ### ECCC_CH
  if(str_detect(name, "T_ECCC_CH_")) {
    species_name <- unlist(str_split(name, "T_ECCC_CH_"))[2]
    ECCC_CH <- c(ECCC_CH, species_name)
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
  if(str_detect(name, "T_IUCN_AMPH_")) {
    species_name <- unlist(str_split(name, "T_IUCN_AMPH_"))[2]
    IUCN_AMPH <- c(IUCN_AMPH, species_name)
  } 
  
  ### IUCN_Birds
  if(str_detect(name, "T_IUCN_BIRD_")) {
    species_name <- unlist(str_split(name, "T_IUCN_BIRD_"))[2]
    IUCN_BIRD <- c(IUCN_BIRD, species_name)
  }
  
  ### IUCN_Mammals
  if(str_detect(name, "T_IUCN_MAMM_")) {
    species_name <- unlist(str_split(name, "T_IUCN_MAMM_"))[2]
    IUCN_MAMM <- c(IUCN_MAMM, species_name)
  }  
  
  ### IUCN_Reptiles
  if(str_detect(name, "T_IUCN_REPT_")) {
    species_name <- unlist(str_split(name, "T_IUCN_REPT_"))[2]
    IUCN_REPT <- c(IUCN_REPT, species_name)
  }  
}


## Create csv ----
## Get the length of the longest vector 
max_ln <- max(c(length(ECCC_SAR), length(ECCC_CH), length(NSC_END),
              length(NSC_SAR), length(NSC_SPP),
              length(IUCN_AMPH), length(IUCN_BIRD),
              length(IUCN_MAMM), length(IUCN_REPT)))

## Build data.frame 
species <- data.frame(ECC_SAR = c(ECCC_SAR, rep(NA, max_ln - length(ECCC_SAR))),
                      ECC_CH = c(ECCC_CH, rep(NA, max_ln - length(ECCC_CH))),
                      NSC_END = c(NSC_END, rep(NA, max_ln - length(NSC_END))),
                      NSC_SAR = c(NSC_SAR, rep(NA, max_ln - length(NSC_SAR))),
                      NSC_SPP = c(NSC_SPP, rep(NA, max_ln - length(NSC_SPP))),
                      IUCN_AMPH = c(IUCN_AMPH, rep(NA, max_ln - length(IUCN_AMPH))),
                      IUCN_BIRD = c(IUCN_BIRD, rep(NA, max_ln - length(IUCN_BIRD))),
                      IUCN_MAMM = c(IUCN_MAMM, rep(NA, max_ln - length(IUCN_MAMM))),
                      IUCN_REPT = c(IUCN_REPT, rep(NA, max_ln - length(IUCN_REPT))))

## Write species data.frame to disk 
write.csv(species, 
          paste0(root_folder, "/Variables/Themes/SPECIES.csv"), 
          row.names = FALSE)

# 6.0 Clear R environment ------------------------------------------------------ 

## End timer
end_time <- Sys.time()
end_time - start_time

rm(list=ls())
gc()
