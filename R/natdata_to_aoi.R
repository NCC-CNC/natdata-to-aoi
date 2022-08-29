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

# Read-in area of interest .tiff (aoi) ----
aoi_path <- "data/aoi/prairie_grasslands.tif" # <------ CHANGE HERE FOR NEW AOI
aoi_1km <- raster(aoi_path) 
aoi_name <- names(aoi_1km)

# Read-in national 1km grid (all of Canada) ----
ncc_1km <- raster("data/national/boundary.tif")
ncc_1km_idx <- ncc_1km
ncc_1km_idx[] <- 1:ncell(ncc_1km_idx) #could restrict to non-NA's is desired

ncc_pu <- ncc_1km[]

# Align aoi to same extent and same number of rows/cols as national grid ----
gdalUtils::align_rasters(unaligned = aoi_path,
                         reference = "data/national/boundary.tif",
                         dstfile = paste0("data/aoi/align/", aoi_name, ".tif"))

# Get numeric vector of aoi planning units 
aoi_pu <- raster("data/aoi/align/R1km_AOI_Aligned.tif")
# aoi_pu <- raster(paste0("data/aoi/align/", aoi_name, ".tif"))#%>%
#   getValues()  %>% .[!is.na(ncc_pu)] 
# # Reclass na to 0
# aoi_pu <- ifelse(!is.na(aoi_pu), 1, 0)
aoi_rij <- prioritizr::rij_matrix(ncc_1km, stack(aoi_pu, ncc_1km_idx)) 
rownames(aoi_rij) <- c("AOI", "Idx")


# Read-in national data ----
SAR <- readRDS("data/national/rij_SAR.rds") 
NSC_END <- readRDS( "data/national/rij_NSC_END.rds")
NSC_SAR <- readRDS( "data/national/rij_NSC_SAR.rds")
NSC_SPP <- readRDS( "data/national/rij_NSC_SPP.rds")

# Matrix multiplication ----
## ECCC Species at risk ----
SAR_rb <- rbind(SAR, aoi_rij)
SAR_rb_red <- SAR_rb[,SAR_rb["AOI",] >0]
SAR_rb_red2 <- SAR_rb_red[rowSums(SAR_rb_red) > 0,]

tt2 <- ncc_1km_idx

rst <- list()

for(ii in 1:(nrow(SAR_rb_red2)-2)){
  tt2[] <- NA
  tt2[SAR_rb_red2["Idx",]] <- SAR_rb_red2[ii,]
  rst[[ii]] <- tt2
}



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
max_ln <- max(c(length(SAR_aoi_pu ), length(NSC_END_aoi_pu )),
              (length(NSC_END_aoi_pu)), length(NSC_SPP_aoi_pu))

## Build data.frame ----
species <- data.frame(ECC_SAR = c(SAR_aoi_pu,rep(NA, max_ln - length(SAR_aoi_pu))),
                      NSC_END = c(NSC_END_aoi_pu,rep(NA, max_ln - length(NSC_END_aoi_pu))),
                      NSC_SAR = c(NSC_SAR_aoi_pu,rep(NA, max_ln - length(NSC_SAR_aoi_pu))),
                      NSC_SPP = c(NSC_SPP_aoi_pu,rep(NA, max_ln - length(NSC_SPP_aoi_pu))))

## Write species data.frame to disk ----
write.csv(species, 
          paste0("data/csv/", names(aoi_1km),"_", species.csv), 
          row.names = FALSE) 

