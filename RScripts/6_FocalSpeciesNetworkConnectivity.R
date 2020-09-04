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
library(sf)
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
focalArea <- st_read(file.path(dataDir, "polygon_projected.shp"))
habitatRaster <- raster(file.path(dataDir, paste0("BLBR_HabitatPatch_buffer", polygonBufferWidth, "km_buffered.tif")))
resistanceRaster <- raster(file.path(dataDir, paste0("BLBR_Resistance_buffer", polygonBufferWidth, "km_buffered.tif")))

# Extract mpg -------------------------------------------------------------
mpg <- MPG(cost = resistanceRaster, patch = habitatRaster)

# Betweenness ------------------------------------------------------------
btwn <- data.frame(patchId = as.numeric(vertex_attr(mpg$mpg, 'patchId')), btwn = betweenness(mpg$mpg, weights = edge_attr(mpg$mpg, 'lcpPerimWeight'), directed = FALSE))
#make a look-up table between patch id and betweenness value
btwn_lookup <- cbind(patchId = btwn$patchId, btwn = 1 / (max(btwn$btwn) - min(btwn$btwn)) * (btwn$btwn - min(btwn$btwn)))
#replace patch ids with betweeness values
btwnRaster <- reclassify(mpg$patchId, btwn_lookup)


# Repeat analysis in focal polygon
# Crop
resistanceRasterFocal <- resistanceRaster %>%
  crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=focalArea) %>% # Clip to focal area
  trim(.) # Trim extra white spaces

habitatRasterFocal <- habitatRaster %>%
  crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  mask(., mask=focalArea) %>% # Clip to focal area
  trim(.) # Trim extra white spaces

mpgFocal <- MPG(cost = resistanceRasterFocal, patch = habitatRasterFocal)

btwnFocal <- data.frame(patchId = as.numeric(vertex_attr(mpgFocal$mpg, 'patchId')), btwn = betweenness(mpgFocal$mpg, weights = edge_attr(mpgFocal$mpg, 'lcpPerimWeight'), directed = FALSE))
#make a look-up table between patch id and betweenness value
btwn_lookup <- cbind(patchId = btwnFocal$patchId, btwn = 1 / (max(btwnFocal$btwn) - min(btwnFocal$btwn)) * (btwnFocal$btwn - min(btwnFocal$btwn)))
#replace patch ids with betweeness values
btwnRasterFocal <- reclassify(mpgFocal$patchId, btwn_lookup)

# Calculate overall EC --------------------------------------------
dist_mat <- shortest.paths(mpgFocal$mpg, weights = edge_attr(mpgFocal$mpg, 'lcpPerimWeight'))

#Natal
#matrix of dispersal kernel
kernel_mat <- exp(-dist_mat / 459)
rm(dist_mat)

#matrix of product node attributes
attribute_mat1 <- as.vector(as.numeric(vertex_attr(mpgFocal$mpg, 'patchArea')))
attribute_mat <- attribute_mat1 %*% t(attribute_mat1)

#matrix of PCvalues
PC_mat <- kernel_mat * attribute_mat
rm(attribute_mat)

PCnum <- sum(PC_mat)
ECNatal <- sqrt(PCnum)


# Save outputs ---------------------------------------------------------
#geotif
writeRaster(btwnRaster, file.path(outDir, paste0("BLBR_Betweenness_", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)
writeRaster(btwnRasterFocal, file.path(outDir, paste0("BLBR_Betweenness_focalArea.tif")), overwrite=TRUE)
writeRaster(mpg$mpgPlot, file.path(outDir, paste0("BLBR_MPG_", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)
writeRaster(mpgFocal$mpgPlot, file.path(outDir, paste0("BLBR_MPG_focalArea.tif")), overwrite=TRUE)
writeRaster(mpg$lcpPerimWeight, file.path(outDir, paste0("BLBR_LinkLength_", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)
writeRaster(mpgFocal$lcpPerimWeight, file.path(outDir, paste0("BLBR_LinkLength_focalArea.tif")), overwrite=TRUE)

