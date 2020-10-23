##################################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping       
# Calculate patch importance metrics for BLBR with corridors and bridges included in the 
# resistance layer
# 09-2020                                    					
#  
#	1. Uses Makurhini R package          
# 		- Function: MK_dPCIIC()  
# 	2. Calculate patch importance metrics for EMBL species    
#    - Using bridge and culvert + resistance layer                                                                                              #     
#  Inputs - 
#     - EMBL species resistance layers (focal)
#     	- including bridge and culvert resistance refined layers (focal)
#     - habitat patches
#	  -dispersal distance 
#
#  Outputs - 
#     - -patch importance and fractions for EMBL with bridges and culverts included
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
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
outDir <- "Results"

  # Input parameters
source(file.path(rawDataDir, "a233_InputParameters.R")) # project level parameters
	polygonBufferWidth
	suitabilityThreshold

## Load data for EMBL species ---------------------------------------------------------
species <- "EMBL"

## Load data
  # Habitat patches
  # Focal 
habitatPatchesFocal <- raster(
						file.path(procDataDir, 
						paste0(species, "_HabitatPatchID_FocalArea.tif")))
resistanceFocal <- raster(
						file.path(procDataDir, 
						paste0(species, "_Resistance_FocalArea.tif")))

  # Bridges & Culverts resistance
resistanceBridgeFocal <- raster(
							file.path(procDataDir,  "GenericBridgeResistanceFocal.tif"))
resistanceCulvertsFocal <- raster(
						file.path(procDataDir,  "GenericCulvertResistanceFocal.tif"))

  # Tabular data
dispersalDistance <- read_csv(
						file.path(
						paste0(rawDataDir, "/Focal Species"), 
						"FocalSpeciesDispersalDistance.csv"))

  #input parameters
maxdist <- suppressWarnings(dispersalDistance[[which(dispersalDistance$Species==species), "Upper", ]])

## Rescale bridge and culvert rasters to have values of 2 for focal species resistance
resistanceBridgeFocalreclass <- reclassify(resistanceBridgeFocal, rcl=cbind(10, 2), include.lowest=T)
resistanceCulvertsFocalreclass <- reclassify(resistanceCulvertsFocal, rcl=cbind(10, 2), include.lowest=T)

## Combine resistance layer with culverts, bridges resistance layer ------------------------------

  #crop to focal extent
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

## Run connectivity metrics using updated resistance values ------------------------------ 
 
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

## Save rasters
  # Focal area 
  #patch importance metrics
  #geotif
writeRaster(PCfocal, 
				file.path(outDir, 
				paste0(species, "_PC_BridgesCulverts_FocalArea.tif")), 
				overwrite=TRUE)
  #geotif
writeRaster(combinedRaster, 
				file.path(procDataDir, 
				paste0(species, "_ResistanceBridgesCulverts_FocalArea.tif")), 
				overwrite=TRUE)


##End script