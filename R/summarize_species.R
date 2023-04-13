# Authors: Dan Wismer & Marc Edwards
#
# Date: April 13th, 2023
#
# Description: Summarizes national species data that intersect your aoi. To be
#              ran after natdata_to_aoi.R
#
# Inputs:  1. root folder (where your "clipped" national layers are located)
#          2. path for species .csv 
#
# Outputs: 1. a csv that lists the species that intersect the aoi 
#
#===============================================================================

# Start timer
start_time <- Sys.time()

# Load packages
library(dplyr)
library(readr)
library(readxl)
library(stringr)
source("R/fct_sci_to_common.R") # <--- CHANGE IF NOT RUNNING SCRIPT FROM GITHUB FOLDER

# INPUT 1: Set root folder ----
root_folder <- "data/output" 

# INPUT 2: Set csv path ----
output_csv <- "/Variables/Themes/SPECIES.csv" 

# Read-in look up tables ----
table_path <- file.path(root_folder, "Variables", "_Tables")
ECCC_SAR_LU <- read_csv(file.path(table_path, "ECCC_SAR_Metadata.csv"))
ECCC_CH_LU <- read_excel(file.path(table_path,  "ECCC_CH_Metadata.xlsx"))
IUCN_LU <- read_csv(file.path(table_path, "IUCN_Metadata.csv"))
NSC_END_LU <- read_excel(file.path(table_path,  "NSC_END_Metadata.xlsx"))
NSC_SAR_LU <- read_excel(file.path(table_path, "NSC_SAR_Metadata.xlsx"))
NSC_SPP_LU <- read_excel(file.path(table_path, "NSC_SPP_Metadata.xlsx"))

# List files in root folder
tiffs <- list.files(
  root_folder, pattern='.tif$', full.names = TRUE, recursive = TRUE
)

# Build empty data.frame (template for metadata.csv)
df <- data.frame(
  Source = character(),
  File = character(),
  Sci = character(),
  Common = character(),
  Area_Km2 = numeric()
)

# Vector of prefixes
prefixes <- c(
  "T_ECCC_CH_", "T_ECCC_SAR_", 
  "T_IUCN_AMPH_", "T_IUCN_BIRD_", "T_IUCN_MAMM_", "T_IUCN_REPT_",
  "T_NSC_END_", "T_NSC_SAR_", "T_NSC_SPP_"
)

# Vector of sources
sources <- c(
  "ECCC CH", "ECCC SAR", 
  "IUCN AMPH", "IUCN BIRD", "IUCN MAMM", "IUCN REPT",
  "NSC END", "NSC SAR", "NSC SPP"
)

# Populate df ----
for (tiff in tiffs) {
  
  ## Get file name ----
  file <- tools::file_path_sans_ext(basename(tiff))
  
  ## Only sort through species Themes
  if (any(startsWith(file, prefixes))) {
    
    ### Get prefix 
    prefix <- prefixes[startsWith(file, prefixes)]
    ### split to get suffix
    suffix <- unlist(str_split(file, prefix))[2] 
    
    ### Get common name ----
    #### ECCC_CH
    if (prefix %in% c("T_ECCC_CH_")) {
      cosewicid <- unlist(str_split(file, "T_ECCC_CH_COSEWICID_"))[2] 
      com_name <- ch_cosewicid_to_name(ECCC_CH_LU, cosewicid, "common")
      #### ECCC_SAR
    } else if (prefix %in% c("T_ECCC_SAR_")) {
      cosewicid <- unlist(str_split(file, "T_ECCC_SAR_COSEWICID_"))[2] 
      com_name <- sar_cosewicid_to_name(ECCC_SAR_LU, cosewicid, "common")
      #### IUCN
    } else if (prefix %in% c("T_IUCN_AMPH_", "T_IUCN_BIRD_", "T_IUCN_MAMM_", "T_IUCN_REPT_")) {
      com_name <- iucn_to_name(IUCN_LU, paste0(suffix, ".tif"))
      #### NSC_END
    } else if (prefix %in% c("T_NSC_END_")) {
      suffix_no_ <- gsub("_", " ", suffix) 
      com_name <- nsc_end_to_name(NSC_END_LU, suffix_no_)
      #### NSC_SAR
    } else if (prefix %in% c("T_NSC_SAR_")) {
      suffix_no_ <- gsub("_", " ", suffix) 
      com_name <- nsc_sar_to_name(NSC_SAR_LU, suffix_no_)
      #### NSC_SPP
    } else if (prefix %in% c("T_NSC_SPP_")) {
      suffix_no_ <- gsub("_", " ", suffix) 
      com_name <- nsc_spp_to_name(NSC_SPP_LU, suffix_no_)
    }
    
    ### Get source ----
    source <- sources[startsWith(file, prefixes)]
    
    ### get sci name ----
    if (prefix %in% c("T_ECCC_CH_")) {
      sci_name <- ch_cosewicid_to_name(ECCC_CH_LU, cosewicid, "sci")
    } else if (prefix %in% c("T_ECCC_SAR_")) {
      sci_name <- sar_cosewicid_to_name(ECCC_SAR_LU, cosewicid, "sci")
    } else {
      sci_name <- gsub("_", " ", suffix) 
    }
    
    ### Get area ----
    species <- raster::raster(tiff)
    area <- raster::cellStats(species, stat = "sum")
    
    ## Append row to data.frame ----
    new_row <- c(source, paste0(file, ".tif"), sci_name, com_name, area)
    df <- structure(rbind(df, new_row), .Names = names(df)) 
  }
}

# OUTPUT 1: write species data.frame to disk 
write.csv(df, file.path(root_folder, output_csv), row.names = FALSE)

# End timer
end_time <- Sys.time()
end_time - start_time
