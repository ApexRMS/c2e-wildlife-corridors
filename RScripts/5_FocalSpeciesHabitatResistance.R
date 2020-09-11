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
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A233 - Cootes to Escarpment"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")

# Input parameters
species <- "EMBL"
polygonBufferWidth <- 20 # In km
suitabilityThreshold <- 60

# Read in data
# Combined LULC layers
LULC <- raster(file.path(outDir, paste0("LULC_", polygonBufferWidth, "km.tif")))
LULC_buffer <- raster(file.path(outDir, paste0("LULC_", polygonBufferWidth, "km_buffered.tif")))
# Focal area polygon
focalArea <- st_read(file.path(outDir, "FocalArea.shp"))
# Study area polygon
studyArea <- st_read(file.path(outDir, "studyarea_20km.shp"))

# Tabular data
crosswalkHabSuit <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "FocalSpeciesHabitatSuitabilityCrosswalk.csv"))
crosswalkResist <- read_csv(file.path(paste0(dataDir, "/Resistance"), "FocalSpeciesResistanceCrosswalk.csv"))
minPatchSize <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "FocalSpeciesMinPatchSize.csv"))
dispersalDistance <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "FocalSpeciesDispersalDistance.csv"))


# Create habitat suitability layer  ---------------------------------------------------------
# Reclassify
suitabilityRaster <- LULC_buffer %>%
  reclassify(., rcl=crosswalkHabSuit[, c("LULC_ID", species)])

# Create habitat patch layer ------------------------------------------------------
patchSizeThreshold <- minPatchSize$MinPatchSizeHa[minPatchSize$Species==species]

#Habitat patch map
habitatRaster <- Which(suitabilityRaster >= suitabilityThreshold) #habitats above suitability threshold
#convert from hectares to m
conversionFromHa <- res(habitatRaster)[1]*res(habitatRaster)[2]*(1/10000)
habitatClump <- clump(habitatRaster)
habitatClumpID <- data.frame(freq(habitatClump))
#Remove clump observations with frequency smaller than minimum habitat patch size (ha)
habitatClumpID <- habitatClumpID[habitatClumpID$count < patchSizeThreshold/conversionFromHa,]
habitatRaster[Which(habitatClump %in% habitatClumpID$value)] <- 0
habitatRasterCont <- clump(habitatRaster) #creates raster with ID's for all patches > min threshold & > suitability threshold


# Create resistance layer ---------------------------------------------------------
# Reclassify
resistanceRasterReclass <- LULC_buffer %>%
  reclassify(., rcl=crosswalkResist[, c("LULC_ID", species)])

#Overlay habitat patches
resistanceRaster <- overlay(resistanceRasterReclass, habitatRaster, fun = function(x,y){return(ifelse(y==1, 1, x))})
#Reclass to assign habitat patches a resistance value = 1 (note that both overlaid values of both 3 and 5 correspond to habitat patches)
#resistanceRaster <- reclassify(resistanceRasterOverlay, rcl = matrix(c(3, 1, 5, 1, 7, 1, 9, 1), ncol=2, byrow = T))

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

habitatRasterContFocal <- habitatRasterCont %>%
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

habitatRasterContStudyArea <- habitatRasterCont %>%
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
writeRaster(suitabilityRasterFocal, file.path(outDir, paste0(species, "_HabitatSuitability_FocalArea.tif")), overwrite=TRUE)
writeRaster(habitatRasterFocal, file.path(outDir, paste0(species, "_HabitatPatch_FocalArea.tif")), overwrite=TRUE)
writeRaster(habitatRasterContFocal, file.path(outDir, paste0(species, "_HabitatPatchCont_FocalArea.tif")), overwrite=TRUE)
writeRaster(resistanceRasterFocal, file.path(outDir, paste0(species, "_Resistance_FocalArea.tif")), overwrite=TRUE)

# Study Area Unbuffered
writeRaster(suitabilityRasterStudyArea, file.path(outDir, paste0(species, "_HabitatSuitability_", polygonBufferWidth, "km.tif")), overwrite=TRUE)
writeRaster(habitatRasterStudyArea, file.path(outDir, paste0(species, "_HabitatPatch_", polygonBufferWidth, "km.tif")), overwrite=TRUE)
writeRaster(habitatRasterContStudyArea, file.path(outDir, paste0(species, "_HabitatPatchCont_",  polygonBufferWidth, "km.tif")), overwrite=TRUE)
writeRaster(resistanceRasterStudyArea, file.path(outDir, paste0(species, "_Resistance_", polygonBufferWidth, "km.tif")), overwrite=TRUE)

# Study Area Buffered
writeRaster(habitatRaster, file.path(outDir, paste0(species, "_HabitatPatch_", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)
writeRaster(resistanceRaster, file.path(outDir, paste0(species, "_Resistance_", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)

