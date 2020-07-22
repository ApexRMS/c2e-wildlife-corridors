# Packages
library(tidyverse)
library(sf)
library(raster)

# Directories
projectDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A233_RBGConnectivity/a233"
dataDir <- file.path(projectDir, "Data/Processed")
resultsDir <- file.path(projectDir, "Results")

# Read in data
studyArea <- st_read(file.path(dataDir, paste0("studyarea_", polygonBufferWidth, "km_unbuffered.shp")))
currMap <- raster(file.path(resultsDir,"try_cum_curmap.tif"))

crs(currMap) <- crs(studyArea)
logcurrMap <- log(currMap)

# Crop to study area
currMap_unbuffered <- currMap %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.)

logcurrMap_unbuffered <- logcurrMap %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.)

writeRaster(currMap, file.path(resultsDir,"try_cum_curmap_bufferred.tif"), overwrite=TRUE)
writeRaster(logcurrMap, file.path(resultsDir,"try_cum_curmap_bufferred_log.tif"), overwrite=TRUE)
writeRaster(currMap_unbuffered, file.path(resultsDir,"try_cum_curmap_unbufferred.tif"), overwrite=TRUE)
writeRaster(logcurrMap_unbuffered, file.path(resultsDir,"try_cum_curmap_unbufferred_log.tif"), overwrite=TRUE)
