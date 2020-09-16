#####################################################################
# a233 Test script - Makurhini library  for 3 focal species
# 09-2020                                       					
#                                                                                                          
#   Inputs (for focal species):
#    -habitat patches
#    -resistance layer
#    -dispersal distance 
#   Function: MK_dPCIIC()
#   Outputs:
#    -patch importance and fractions
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

# Workspace ---------------------------------------------------------

## Packages
  #To initially install
#require(devtools)
#install_github("connectscape/Makurhini", dependencies = TRUE, upgrade = "never", force=TRUE)

library(Makurhini)

  # Load Additional packages
library(tidyverse)
library(raster)
library(sf)


 # Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A233 - Cootes to Escarpment"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")


## Run for focal species------
species <- "ODVI"		#BLBR EMBL
polygonBufferWidth <- 20 # In km
suitabilityThreshold <- 60

## Load data
  # Habitat patches
  # Focal 
habitatPatchesFocal <- raster(file.path(outDir, paste0(species, "_HabitatPatchCont_FocalArea.tif")))
resistanceFocal <- raster(file.path(outDir, paste0(species, "_Resistance_FocalArea.tif")))

 # 20km
habitatPatches20km <- raster(file.path(outDir, paste0(species,  "_HabitatPatchCont_", polygonBufferWidth, "km.tif")))
resistance20km <- raster(file.path(outDir, paste0(species,  "_Resistance_", polygonBufferWidth, "km_buffered.tif")))

  # Tabular data
dispersalDistance <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "FocalSpeciesDispersalDistance.csv"))

## Input parameters
maxdist <- dispersalDistance[[which(dispersalDistance$Species==species), "Upper", ]]


## Test of IIC and fractions measurements for focal species-----

## Run with least-cost distances, using resistance values. 
  #Focal area
  # IIC
#IICfocal <- MK_dPCIIC(nodes = habitatPatchesFocal, 
#                distance = list(type = "least-cost", 
#                				resistance = resistanceFocal, 
#                				mask =  habitatPatchesFocal),
#                attribute = NULL,
#                metric = "IIC", 
#                distance_thresholds = maxdist,
#                rasterparallel = T) 


  # PC & fractions
PCfocal <- MK_dPCIIC(nodes = habitatPatchesFocal, 
                distance = list(type = "least-cost", 
                				resistance = resistanceFocal, 
                				mask =  habitatPatchesFocal),
                attribute = NULL,
                metric = "PC", 
                probability = 0.95, 
                distance_thresholds = maxdist,
                rasterparallel = T)
plot(PCfocal[["dPCflux"]]) #Test plot dPCflux


##20 km
  # PC & fractions
PC20 <- MK_dPCIIC(nodes = habitatPatches20km, 
                distance = list(type = "least-cost", 
                				resistance = resistance20km, 
                				mask =  habitatPatches20km),
                attribute = NULL,
                metric = "PC", 
                probability = 0.95, 
                distance_thresholds = maxdist,
                rasterparallel = T)
plot(PCfocal[["dPCflux"]]) #Test plot dPCflux


## Save rasters
  #Focal area
#writeRaster(IICfocal, file.path(outDir, paste0(species, "_IIC_FocalArea.tif")), overwrite=TRUE)
writeRaster(PCfocal, file.path(outDir, paste0(species, "_PC_FocalArea.tif")), overwrite=TRUE)

  # Study Area Unbuffered
#writeRaster(IIC20km, file.path(outDir, paste0(species, "_IIC_", polygonBufferWidth, "km.tif")), overwrite=TRUE)
writeRaster(PC20km, file.path(outDir, paste0(species, "_PC_", polygonBufferWidth, "km.tif")), overwrite=TRUE)
