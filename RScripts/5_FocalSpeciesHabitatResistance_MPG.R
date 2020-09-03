#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping            #
# Create habitat and resistance maps for focal species              #
#                                                                   #
# Inputs:                                                           #
#    - The combined LULC layer produced by 1_StudyAreaLULC.R        #
#    - The LULC/habitat suitability crosswalk                       #
# Outputs:                                                          #
#    - Habitat suitability map                                      #
#    - Habitat patch map                                            #
#    - Resistance map                                               #
#                                                                   #
# Script created by Bronwyn Rayfield for ApexRMS                    #
#####################################################################

# Workspace ---------------------------------------------------------
# Packages
library(tidyverse)
library(raster)
library(sf)

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

# Directories
projectDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A233_RBGConnectivity/a233"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")

# Input parameters
polygonBufferWidth <- 20 # In km

# Read in data
# Combined LULC layers
LULC <- raster(file.path(outDir, paste0("LULC_", polygonBufferWidth, "km_buffered.tif")))
LULC_unbuffered <- raster(file.path(outDir, paste0("LULC_", polygonBufferWidth, "km_unbuffered.tif")))
# Focal area polygon
focalArea <- st_read(file.path(dataDir, "polygon_projected.shp"))
# Study area polygon
studyArea <- st_read(file.path(dataDir, "studyarea_20km_unbuffered.shp"))

# Resistance crosswalk
crosswalkHabSuit <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "HabitatSuitabilityCrosswalk.csv"))
crosswalkResist <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "FocalSpeciesResistanceCrosswalk.csv"))

# Create habitat suitability layer  ---------------------------------------------------------
# Reclassify
suitabilityRaster <- LULC %>%
  reclassify(., rcl=crosswalkHabSuit[, c("LULC_ID", "BLBR")])

# Create habitat patch layer ------------------------------------------------------
suitabilityThreshold <- 60
patchSizeThreshold <- 1

#Habitat patch map
habitatRaster <- Which(suitabilityRaster >= suitabilityThreshold)
#convert from hectares to m
conversionFromHa <- res(habitatRaster)[1]*res(habitatRaster)[2]*(1/10000)
habitatClump <- clump(habitatRaster)
habitatClumpID <- data.frame(freq(habitatClump))
#Remove clump observations with frequency smaller than minimum habitat patch size (ha)
habitatClumpID <- habitatClumpID[habitatClumpID$count < patchSizeThreshold/conversionFromHa,]
habitatRaster[Which(habitatClump %in% habitatClumpID$value)] <- 0

# Create resistance layer ---------------------------------------------------------
# Reclassify
resistanceRasterReclass <- LULC %>%
  reclassify(., rcl=crosswalkResist[, c("LULC_ID", "BLBR")])
#Overlay habitat patches
resistanceRasterOverlay <- overlay(resistanceRasterReclass, habitatRaster, fun = function(x,y){return(x+y)})
#Reclass to assign habitat patches a resistance value = 1 (note that both overlaid values of both 3 and 5 correspond to habitat patches)
resistanceRaster <- reclassify(resistanceRasterOverlay, rcl = matrix(c(3, 1, 5, 1), ncol=2, byrow = T))

# Crop to focal region and study area
# Focal Area
suitabilityRasterFocal <- suitabilityRaster %>%
  crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=focalArea) %>% # Clip to focal area
  trim(.) # Trim extra white spaces  

habitatRasterFocal <- habitatRaster %>%
  crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=focalArea) %>% # Clip to focal area
  trim(.) # Trim extra white spaces

resistanceRasterFocal <- resistanceRaster %>%
  crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=focalArea) %>% # Clip to focal area
  trim(.) # Trim extra white spaces

# Study Area
suitabilityRasterStudyArea <- suitabilityRaster %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=studyArea) %>% # Clip to focal area
  trim(.) # Trim extra white spaces  

habitatRasterStudyArea <- habitatRaster %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=studyArea) %>% # Clip to focal area
  trim(.) # Trim extra white spaces

resistanceRasterStudyArea <- resistanceRaster %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=studyArea) %>% # Clip to focal area
  trim(.) # Trim extra white spaces


# Save outputs ---------------------------------------------------------
#geotif
# Focal Area
writeRaster(suitabilityRasterFocal, file.path(outDir, paste0("BLBR_HabitatSuitability_FocalArea.tif")), overwrite=TRUE)
writeRaster(habitatRasterFocal, file.path(outDir, paste0("BLBR_HabitatPatch_FocalArea.tif")), overwrite=TRUE)
writeRaster(resistanceRasterFocal, file.path(outDir, paste0("BLBR_Resistance_FocalArea.tif")), overwrite=TRUE)
# Study Area Unbuffered
writeRaster(suitabilityRasterStudyArea, file.path(outDir, paste0("BLBR_HabitatSuitability_buffer", polygonBufferWidth, "km.tif")), overwrite=TRUE)
writeRaster(habitatRasterStudyArea, file.path(outDir, paste0("BLBR_HabitatPatch_buffer", polygonBufferWidth, "km.tif")), overwrite=TRUE)
writeRaster(resistanceRasterStudyArea, file.path(outDir, paste0("BLBR_Resistance_buffer", polygonBufferWidth, "km.tif")), overwrite=TRUE)
# Study Area Buffered
writeRaster(suitabilityRaster, file.path(outDir, paste0("BLBR_HabitatSuitability_buffer", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)
writeRaster(habitatRaster, file.path(outDir, paste0("BLBR_HabitatPatch_buffer", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)
writeRaster(resistanceRaster, file.path(outDir, paste0("BLBR_Resistance_buffer", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)

