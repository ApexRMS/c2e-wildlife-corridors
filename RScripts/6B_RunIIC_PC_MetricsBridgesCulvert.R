##################################################################################
# a233 Cootes to Escarpment
# 09-2020                                       					
#  Calculate connectivity metrics for test species with Bridges, Culverts
#  
#  1. Calculate PC metrics for EMBL species    
#    - Using bridge and culvert + resistance layer                                                                                              #     
#  Inputs - 
#     - focal species resistance layers (focal)
#     - bridge and culvert resistance refined layers (focal)
#     - habitat patches w IDs
#
#  Outputs - 
#     - PC raster layer for EMBL with bridges and culverts included
#                                                                   
# Script by C Tucker for ApexRMS 									
##################################################################################

## Workspace ---------------------------------------------------------

  # Load packages
library(Makurhini)  
library(tidyverse)
library(raster)
library(sf)

  # Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A233 - Cootes to Escarpment"
rawDataDir <- file.path(projectDir, "Data/Raw")
procDataDir <- file.path(projectDir, "Data/Processed")
outDir <- file.path(projectDir, "Results")

## Run for EMBL species------
species <- "EMBL"
polygonBufferWidth <- 20 # In km
suitabilityThreshold <- 60

## Load data
  # Habitat patches
  # Focal 
habitatPatchesFocal <- raster(file.path(procDataDir, paste0(species, "_HabitatPatchID_FocalArea.tif")))
resistanceFocal <- raster(file.path(procDataDir, paste0(species, "_Resistance_FocalArea.tif")))

  # + Bridges & Culverts resistance, focal scale
resistanceBridgeFocal <- raster(file.path(procDataDir,  "GenericBridgeResistanceFocal.tif"))
resistanceCulvertsFocal <- raster(file.path(procDataDir,  "GenericCulvertResistanceFocal.tif"))

  # Tabular data
dispersalDistance <- read_csv(file.path(paste0(rawDataDir, "/Focal Species"), "FocalSpeciesDispersalDistance.csv"))

  #input parameters
maxdist <- suppressWarnings(dispersalDistance[[which(dispersalDistance$Species==species), "Upper", ]])


## Rescale bridge and culvert rasters to have values of 2 for focal species resistance
resistanceBridgeFocalreclass <- reclassify(resistanceBridgeFocal, rcl=cbind(10, 2), include.lowest=T)
resistanceCulvertsFocalreclass <- reclassify(resistanceCulvertsFocal, rcl=cbind(10, 2), include.lowest=T)

## Combine resistance layer with culverts, bridges resistance layer
  #crop to same extent
resistanceRasterFocal <- resistanceFocal %>%
  crop(., extent(habitatPatchesFocal), snap="out") # Crop to focal area extent
 
bridgeRasterFocal <- resistanceBridgeFocalreclass %>%
  crop(., extent(habitatPatchesFocal), snap="out") 

culvertsRasterFocal <- resistanceCulvertsFocalreclass %>%
  crop(., extent(habitatPatchesFocal), snap="out") 
 

  # Combine rasters by selecting min resistance value per cell
combinedRaster <- resistanceRasterFocal %>%
					stack(., bridgeRasterFocal, culvertsRasterFocal) %>%
					min(., na.rm=TRUE)

par(mfrow=c(1, 2))
plot(resistanceRasterFocal)
plot(combinedRaster)

## Run connectivity metrics using updated resistance values. 
 
  # PC & fractions
PCfocal <- MK_dPCIIC(nodes = habitatPatchesFocal, 
                distance = list(type = "least-cost", 
                				resistance = combinedRaster, 
                				mask =  habitatPatchesFocal),
                attribute = NULL,
                metric = "PC", 
                probability = 0.95, 
                distance_thresholds = maxdist,
                rasterparallel = T)
# plot(PCfocal[["dPCflux"]]) #Test plot dPCflux


## Save rasters
  # Focal area
writeRaster(PCfocal, file.path(outDir, paste0(species, "_PC_FocalAreaCulvertsBridges.tif")), overwrite=TRUE)
writeRaster(combinedRaster, file.path(procDataDir, paste0(species, "_RefinedResistanceFocalArea.tif")), overwrite=TRUE)

