#####################################################################
# a233 Test script - Makurhini library 
# 09-2020                                       					#
#                                                                   #                                       #
#                                                                   #
# Script by C Tucker for ApexRMS 									#
#####################################################################

# Workspace ---------------------------------------------------------

## Packages
  #To initially install
#require(devtools)
#install_github("connectscape/Makurhini", dependencies = TRUE, upgrade = "never", force=TRUE)

library(Makurhini)

## Test of IIC and fractions measurements -----

  #Load data
data("vegetation_patches", package = "Makurhini")
nrow(vegetation_patches)
# 142

  #IIC
IIC <- MK_dPCIIC(nodes = vegetation_patches, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "IIC", distance_thresholds = 10000) #10 km

plot(IIC["dIIC"], breaks="jenks")


  # PC & fractions
PC <- MK_dPCIIC(nodes = vegetation_patches, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "PC", probability = 0.05,
                distance_thresholds = 10000)  
                
plot(PC["dPCflux"], breaks="jenks")



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


## Run for EMBL, Blandings turtle------
species <- as.character("EMBL")
polygonBufferWidth <- 20 # In km
suitabilityThreshold <- 60

  # Load data
  # Habitat patches
  # Focal 
habitatPatchesFocal <- raster(file.path(outDir, paste0(species, "_HabitatPatchCont_FocalArea.tif")))
resistanceFocal <- raster(file.path(outDir, paste0(species, "_Resistance_FocalArea.tif")))
habitatSuitabilityFocal <- raster(file.path(outDir, paste0(species, "_HabitatSuitability_FocalArea.tif")))

 # 20km
habitatPatches20km <- raster(file.path(outDir, paste0(species,  "_HabitatPatchCont_", polygonBufferWidth, "km.tif")))
resistance20km <- raster(file.path(outDir, paste0(species,  "_Resistance_", polygonBufferWidth, "km_buffered.tif")))
habitatSuitability20km <- raster(file.path(outDir, paste0(species,  "_HabitatSuitability_", polygonBufferWidth, "km.tif")))


# Tabular data
crosswalkHabSuit <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "FocalSpeciesHabitatSuitabilityCrosswalk.csv"))
crosswalkResist <- read_csv(file.path(paste0(dataDir, "/Resistance"), "FocalSpeciesResistanceCrosswalk.csv"))
minPatchSize <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "FocalSpeciesMinPatchSize.csv"))
dispersalDistance <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "FocalSpeciesDispersalDistance.csv"))
#plot(habitatPatchesFocal)

  #Attributes
maxdist <- dispersalDistance[dispersalDistance$Species==species, "Upper"]

# Extract habitat suitability for habitat patches
habitatSuitabilityFocalCrop <- mask(crop(habitatSuitabilityFocal, habitatPatchesFocal), habitatPatchesFocal)


  #IIC
IIC <- MK_dPCIIC(nodes = habitatPatchesFocal, 
				area_unit = "m2",
				attribute = NULL,
                distance = list(type= "least-cost", resistance = resistanceFocal),
                metric = "IIC", 
                probability = 0.95, 
                distance_thresholds = maxdist) #based on max distance 
                plot(IIC)

plot(IIC[["dIIC"]]) #Test plot dIIC


  # PC & fractions
PC <- MK_dPCIIC(nodes = habitatPatchesFocal, attribute = NULL,
                distance = list(type= "least-cost", resistance = resistanceFocal),
                metric = "PC", probability = 0.95,
                distance_thresholds = maxdist)  
plot(PC[["dPCflux"]]) #Test plot dPCflux

#class      : RasterStack 
#dimensions : 750, 582, 436500, 5  (nrow, ncol, ncell, nlayers)
#resolution : 15, 15  (x, y)
#extent     : 1340475, 1349205, 11863365, 11874615  (xmin, xmax, ymin, ymax)
#crs        : +proj=lcc +lat_1=44.5 +lat_2=53.5 +lat_0=0 +lon_0=-85 +x_0=930000 +y_0=6430000 +ellps=GRS80 +units=m +no_defs 
#names      :          Id,         dPC,    dPCintra,     dPCflux, dPCconnector 
#min values : 1.000000000, 1.231267745, 0.006023568, 1.225244177,  0.000000000 
#max values :    32.00000,    85.82702,    42.20662,    43.62041,      0.00000 



## Run for BLBR------
species <- as.character("BLBR")

## Load data
BLBR <- raster(file.path(outDir, paste0(species, "_Resistance_FocalArea.tif")))


## Run for BLBR
  #IIC
IIC <- MK_dPCIIC(nodes = BLBR, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "IIC", distance_thresholds = 10000) #10 km
                plot(IIC)
plot(IIC[["dIIC"]]) #Test plot dIIC
                
#class      : RasterStack 
#dimensions : 750, 582, 436500, 5  (nrow, ncol, ncell, nlayers)
#resolution : 15, 15  (x, y)
#extent     : 1340475, 1349205, 11863365, 11874615  (xmin, xmax, ymin, ymax)
#crs        : +proj=lcc +lat_1=44.5 +lat_2=53.5 +lat_0=0 +lon_0=-85 +x_0=930000 +y_0=6430000 +ellps=GRS80 +units=m +no_defs 
#names      :           Id,         dIIC,    dIICintra,     dIICflux, dIICconnector 
#min values :  1.000000000,  0.942409948,  0.005587353,  0.936822595,   0.000000000 
#max values : 3.200000e+01, 5.411643e+01, 1.842403e+01, 3.569240e+01,  7.105427e-15 



  # PC & fractions
PC <- MK_dPCIIC(nodes = BLBR, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "PC", probability = 0.05,
                distance_thresholds = 10000)  
plot(PC[["dPCflux"]]) #Test plot dPCflux

#class      : RasterStack 
#dimensions : 750, 582, 436500, 5  (nrow, ncol, ncell, nlayers)
#resolution : 15, 15  (x, y)
#extent     : 1340475, 1349205, 11863365, 11874615  (xmin, xmax, ymin, ymax)
#crs        : +proj=lcc +lat_1=44.5 +lat_2=53.5 +lat_0=0 +lon_0=-85 +x_0=930000 +y_0=6430000 +ellps=GRS80 +units=m +no_defs 
#names      :           Id,          dPC,     dPCintra,      dPCflux, dPCconnector 
#min values :  1.000000000,  0.754015632,  0.004329759,  0.749685874,  0.000000000 
#max values : 3.200000e+01, 5.834430e+01, 1.427717e+01, 4.406713e+01, 7.660539e-15 


