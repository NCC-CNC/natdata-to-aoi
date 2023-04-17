# Extract national raster layers to AOI

The `natdata_to_aoi.R` script extracts national data (1km cell size) to an area of interest (aoi) raster.

:warning: **Data needed to run this script is not packaged in this repo.**
 
## Inputs:
1. Root folder output path
2. Input path to your AOI.tif

## Outputs:
1. National raster layer for each variable that intersects the aoi
2. Csv that lists the species that intersect the aoi

## National Data:
### Themes ----

* Environment and Climate Change Canada Critical Habitat (ECCC_CH)
* Environment and Climate Change Canada Species at Risk (ECCC_SAR)
* International Union for Conservation of Nature Amphibians (IUCN_AMPH)
* International Union for Conservation of Nature Birds (IUCN_BIRD)
* International Union for Conservation of Nature Mammals (IUCN_MAMM)
* International Union for Conservation of Nature Reptiles (IUCN_REPT)
* Nature Serve Canada Species at Risk (NSC_SAR)
* Nature Serve Canada Endemics (NSC_END)
* Nature Serve Canada Common Species (NSC_SPP)
* Forest Land Cover
* Forest Land Use
* Wetlands
* Grasslands
* Lakes
* Rivers
* Shoreline

### Weights ----

* Carbon storage
* Carbon potential
* Climate forward velocity
* Climate refugia
* Climate extremes
* Connectivity
* Human Footprint Index
* Key Biodiversity Areas
* Recreation
* Freshwater

### Includes ----

* Existing Conservation (CPCAD + NCC Fee Simple and Conservation Agreements)

