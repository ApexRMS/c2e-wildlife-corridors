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

# Workspace ---------------------------------------------------------
# Packages
library(tidyverse)
library(raster)

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

# Directories
projectDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A233_RBGConnectivity/a233"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")

# Input parameters
polygonBufferWidth <- 20 # In km

# Read in data
# Resistance layer
resistance <- raster(file.path(outDir, paste0("Resistance_buffer", polygonBufferWidth, "km_buffered.tif")))

# Create binary study area raster
studyArea <- resistance
studyArea[studyArea>0] <- 1

boundary <- boundaries(studyArea, type='inner', directions=8, asNA=T)

focalNodes5 <- sampleRandom(boundary, size=5, asRaster=T)
focalNodesID5 <- clump(focalNodes5)


# y[1,2448]<-1
# y[nrow(y),1990]<-2
# y[1,]<-1
# y[nrow(y),]<-2


y[y==0]<-NA
writeRaster(focalNodesID5, file.path(outDir, paste0("FocalNode_try5_", polygonBufferWidth, "km_buffered.asc")), overwrite=T)
writeRaster(y, file.path(outDir, paste0("FocalNode_test_", polygonBufferWidth, "km_buffered.tif")), overwrite=T)
writeRaster(y, file.path(outDir, paste0("FocalNode_test_", polygonBufferWidth, "km_buffered.asc")), overwrite=T)

unique(y)