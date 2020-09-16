##################################################################################
# a233 Run Makurhini library to calculate connectivity metrics for test species with Bridges, Culverts
# 09-2020                                       					
#  1. Calculate IIC and PC metrics for EMBL species                                                                                                      
#  2. Calculate IIC and PC metrics for EMBL species  
#     - using bridge and culvert + resistance layer 
#
#  Inputs - 
#     - generic resistance layers (focal)
#     - generic resistance refined layers (focal)
#     - habitat patches w IDs
#
#  Outputs - 
#     - IIC and PC raster layers
#                                                                   
# Script by C Tucker for ApexRMS 									
##################################################################################

# Workspace ---------------------------------------------------------

## Packages
  #To initially install
#require(devtools)
#install_github("connectscape/Makurhini", dependencies = TRUE, upgrade = "never", force=TRUE)

library(Makurhini)


##--------------------------------------------------------------------------
## Test of IIC and fractions measurements for EMBL focal species-----

  # Load Additional packages
library(tidyverse)
library(raster)
library(sf)

 # Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A233 - Cootes to Escarpment"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")


## Run for focal species------
species <- "EMBL"
polygonBufferWidth <- 20 # In km
suitabilityThreshold <- 60

## Load data
  # Habitat patches
  # Focal 
habitatPatchesFocal <- raster(file.path(outDir, paste0(species, "_HabitatPatchCont_FocalArea.tif")))
resistanceFocal <- raster(file.path(outDir, paste0(species, "_Resistance_FocalArea.tif")))

  # + Bridges & Culverts resistance, focal scale
resistanceBridgeFocal <- raster(file.path(outDir,  "GenericSpeciesBridgeResistanceFocal.tif"))
resistanceCulvertsFocal <- raster(file.path(outDir,  "GenericSpeciesCulvertResistanceFocal.tif"))

  # Tabular data
dispersalDistance <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "FocalSpeciesDispersalDistance.csv"))

  #input parameters
maxdist <- dispersalDistance[[which(dispersalDistance$Species==species), "Upper", ]]


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
resistanceStack <- stack(resistanceRasterFocal, bridgeRasterFocal, culvertsRasterFocal)
combinedRaster  <- min(resistanceStack, na.rm=TRUE)

par(mfrow=c(1, 2))
plot(resistanceRasterFocal)
plot(combinedRaster)

## Run with least-cost distances, using resistance values. 
  #Focal area

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
plot(PCfocal[["dPCflux"]]) #Test plot dPCflux


## Save rasters
  #Focal area
writeRaster(PCfocal, file.path(outDir, paste0(species, "_PC_FocalAreaCulvertsBridges.tif")), overwrite=TRUE)

writeRaster(combinedRaster, file.path(outDir, paste0(species, "ResistanceFocalAreaCulvertsBridges.tif")), overwrite=TRUE)