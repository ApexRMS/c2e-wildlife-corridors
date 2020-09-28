# Packages
library(tidyverse)
library(sf)
library(raster)

# Directories
projectDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A233_RBGConnectivity/a233"
dataDir <- file.path(projectDir, "Data/Processed")
resultsDir <- file.path(projectDir, "Results")

analysisName <- "Generic_CulvertBridge_5pairwise"

# Read in data
studyArea <- st_read(file.path(dataDir, paste0("StudyArea_20km.shp")))
focalArea <- st_read(file.path(dataDir, "FocalArea.shp"))
currMap_buffer <- raster(file.path(resultsDir,paste0(analysisName, "_cum_curmap.asc")))

crs(currMap_buffer) <- crs(studyArea)
logcurrMap_buffer <- log(currMap_buffer)

# Crop to study area
currMap <- currMap_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.)

logcurrMap <- logcurrMap_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.)

# Crop to focal area
currMap_focalArea <- currMap_buffer %>%
  crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=focalArea) %>% # Clip to focal area
  trim(.)

logcurrMap_focalArea <- logcurrMap_buffer %>%
  crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=focalArea) %>% # Clip to focal area
  trim(.)

# Save geotifs
# Study Area buffer
writeRaster(currMap_buffer, file.path(resultsDir, paste0(analysisName, "_cum_curmap_20km_buffer.tif")), overwrite=TRUE)
writeRaster(logcurrMap_buffer, file.path(resultsDir, paste0(analysisName, "_cum_curmap_log_20km_buffer.tif")), overwrite=TRUE)
# Study Area
writeRaster(currMap, file.path(resultsDir, paste0(analysisName, "_cum_curmap_20km.tif")), overwrite=TRUE)
writeRaster(logcurrMap, file.path(resultsDir, paste0(analysisName, "_cum_curmap_log_20km.tif")), overwrite=TRUE)
# Focal Area
writeRaster(currMap_focalArea, file.path(resultsDir, paste0(analysisName, "_cum_curmap_focalArea.tif")), overwrite=TRUE)
writeRaster(logcurrMap_focalArea, file.path(resultsDir, paste0(analysisName, "_cum_curmap_log_focalArea.tif")), overwrite=TRUE)
