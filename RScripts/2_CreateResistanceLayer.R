#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping            #
# Create resistance surface                                         #
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
      # Combined LULC layers
LULC <- raster(file.path(outDir, paste0("LULC_", polygonBufferWidth, "km_buffered.tif")))
LULC_unbuffered <- raster(file.path(outDir, paste0("LULC_", polygonBufferWidth, "km_unbuffered.tif")))

      # Resistance crosswalk
crosswalk <- read_csv(file.path(paste0(dataDir, "/Resistance"), "ResistanceCrosswalk.csv"))

# Create resistance layer  ---------------------------------------------------------
# Reclassify
resistance <- LULC %>%
  reclassify(., rcl=crosswalk[, c("Destination_ID", "Resistance")])

resistance_unbuffered <- LULC_unbuffered %>%
  reclassify(., rcl=crosswalk[, c("Destination_ID", "Resistance")])

# Save outputs ---------------------------------------------------------
#geotif
writeRaster(resistance, file.path(outDir, paste0("Resistance_buffer", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)
writeRaster(resistance_unbuffered, file.path(outDir, paste0("Resistance_buffer", polygonBufferWidth, "km_unbuffered.tif")), overwrite=TRUE)
#asc
writeRaster(resistance, file.path(outDir, paste0("Resistance_buffer", polygonBufferWidth, "km_buffered.asc")), overwrite=TRUE)
writeRaster(resistance_unbuffered, file.path(outDir, paste0("Resistance_buffer", polygonBufferWidth, "km_unbuffered.asc")), overwrite=TRUE)



