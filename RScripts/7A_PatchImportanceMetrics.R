#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping       
# Calculate patch importance metrics for 3 focal species
# 09-2020                                       					
#                             
#	1. Uses Makurhini R package          
# 		- Function: MK_dPCIIC()   
#                                             
#   Inputs (for focal species):
#    -habitat patches
#    -resistance layer
#    -dispersal distance 
#  
#   Outputs:
#    -patch importance and fractions
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

## Workspace ---------------------------------------------------------

  # Packages
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
	specieslist
	polygonBufferWidth
	suitabilityThreshold

## Calculate patch importance for focal species -----------------------------------------

  # Loop through species, generating Patch Importance data and saving outputs

for (i in specieslist){
	
	species <- i

## Load data for focal species
  # Focal area
habitatPatchesFocal <- raster(
						file.path(procDataDir, 
						paste0(species, "_HabitatPatchID_FocalArea.tif")))
resistanceFocal <- raster(
						file.path(procDataDir, 
						paste0(species, "_Resistance_FocalArea.tif")))

  # Tabular data
dispersalDistance <- read_csv(
						file.path(
						paste0(rawDataDir, "/Focal Species"), 
						"FocalSpeciesDispersalDistance.csv"))

 # Species dispersal maximal distance 
maxdist <- suppressWarnings(dispersalDistance[[which(dispersalDistance$Species == species), "Upper", ]])

## Test of IIC and fractions measurements for focal species ----------------------------------------

  # Focal area
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

## Save output files
  # Focal area
  #geotif
writeRaster(PCfocal, 
				file.path(outDir, 
				paste0(species, "_PC_FocalArea.tif")), 
				overwrite=TRUE)  

rm(PCfocal)

} #end loop

# End script