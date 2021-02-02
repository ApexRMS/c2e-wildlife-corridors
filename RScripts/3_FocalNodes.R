#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping            #
# Create focal nodes                                                #
#                                                                   #
# Inputs:                                                           #
#    - The combined LULC layer produced by 1_StudyAreaLULC.R        #                                                                       #
#    - The LULC/resistance crosswalk                                #
# Outputs:                                                          #
#    - Resistance layer                                             #
#                                                                   #
# Script created by Bronwyn Rayfield and Chloé Debeyser for ApexRMS #
#####################################################################

## Workspace ---------------------------------------------------------
  # Packages
library(tidyverse)
library(raster)

  # Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

  # Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"

## Input parameters
source(file.path("Data", "a233_InputParams.R")) # project level parameters
	polygonBufferWidth


## Read in data---------------------------------------------------------
  # Resistance layer
studyArea <- raster(
				file.path(procDataDir, 
				paste0("Generic_Resistance_", polygonBufferWidth, "km_buffer.tif")))


## Generate nodes ---------------------------------------------------------

  # Create binary study area raster
studyArea[studyArea > 0] <- 1
  
  #ID boundary 
boundary <- boundaries(studyArea, type='inner', directions=8, asNA=T)

focalNodes50 <- sampleRandom(boundary, size=50, asRaster=T)
focalNodesID50 <- clump(focalNodes50)



writeRaster(focalNodesID50, file.path(procDataDir, paste0("FocalNode_50_", polygonBufferWidth, "km_buffered.asc")), overwrite=T)
writeRaster(y, file.path(procDataDir, paste0("FocalNode_test_", polygonBufferWidth, "km_buffered.tif")), overwrite=T)
writeRaster(y, file.path(procDataDir, paste0("FocalNode_test_", polygonBufferWidth, "km_buffered.asc")), overwrite=T)



###??

# y[1,2448]<-1
# y[nrow(y),1990]<-2
# y[1,]<-1
# y[nrow(y),]<-2


y[y==0]<-NA
unique(y)