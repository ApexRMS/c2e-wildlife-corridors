#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping            #
# Create habitat and resistance maps for focal species              #
#                                                                   #
# Inputs:                                                           #
#    - Habitat suitability map                                      #
#    - Resistance map                                               #
#    - Dispersal parameters                                         #
# Outputs:                                                          #
#    - Habitat network map                                          #
#    - Betweenness map                                              #
#                                                                   #
# Script created by Bronwyn Rayfield for ApexRMS                    #
#####################################################################

# Workspace ---------------------------------------------------------
# Packages
library(tidyverse)
library(raster)
library(igraph)
library(grainscape)

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

# Directories
projectDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A233_RBGConnectivity/a233"
dataDir <- file.path(projectDir, "Data/Processed")
outDir <- file.path(projectDir, "Data/Processed")

# Input parameters
polygonBufferWidth <- 20 # In km

# Read in data
# Combined LULC layers
habitatRaster <- raster(file.path(outDir, paste0("BLBR_HabitatPatch_buffer", polygonBufferWidth, "km_buffered.tif")))
resistanceRaster <- raster(file.path(outDir, paste0("BLBR_Resistance_buffer", polygonBufferWidth, "km_buffered.tif")))

# Extract mpg -------------------------------------------------------------
mpg <- MPG(cost = resistanceRaster, patch = habitatRaster)

# Betweenness ------------------------------------------------------------
btwn <- data.frame(patchId = as.numeric(vertex_attr(mpg$mpg, 'patchId')), btwn = betweenness(mpg$mpg, weights = edge_attr(mpg$mpg, 'lcpPerimWeight'), directed = FALSE))
#make a look-up table between patch id and betweenness value
btwn_lookup <- cbind(patchId = btwn$patchId, btwn = 1 / (max(btwn$btwn) - min(btwn$btwn)) * (btwn$btwn - min(btwn$btwn)))
#replace patch ids with betweeness values
btwnRaster <- reclassify(mpg$patchId, btwn_lookup)


