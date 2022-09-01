### Overview
The `natdata_to_aoi.R` script extracts national data (1km cell size) to an area of interest (aoi).

**Note:**
data needed to run this script is not packaged in this repo. 

### Inputs:
1. an output folder path
2. an aoi raster

### Outputs:
1. a raster layer (1km cell size) for each variable that intersects the aoi
2. a csv that lists the species that intersect the aoi

### National Data:
**Species**

* Environment and Climate Change Canada Species at Risk (ECCC_SAR)
* Nature Serve Canada Species at Risk (NSC_SAR)
* Nature Serve Canada Endemics (NSC_END)
* Nature Serve Canada Common Species (NSC_SPP)
* International Union for Conservation of Nature Amphibians (IUCN_Amphibians)
* International Union for Conservation of Nature Birds (IUCN_Birds)
* International Union for Conservation of Nature Mammals (IUCN_Mammals)
* International Union for Conservation of Nature Reptiles (IUCN_Reptiles)

**Conservation Variables**

* Carbon
* Climate
* Forest
* Grassland
* Protected
* Recreation
* Water
* Wetlands
