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

  # Functions
rescaleR <- function(x, new.min = 0, new.max = 1) {
   x.min = suppressWarnings(min(x, na.rm=TRUE))
   x.max = suppressWarnings(max(x, na.rm=TRUE))
   new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

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

 # Species dispersal maximal/or median distance 
maxdist <- suppressWarnings(dispersalDistance[[which(dispersalDistance$Species == species), "Upper", ]])

prob <- ifelse(species == "EMBL", 0.95, 0.5)

## Test of IIC and fractions measurements for focal species ----------------------------------------

  # Focal area
  # PC & fractions
PCfocal <- MK_dPCIIC(nodes = habitatPatchesFocal, 
                	 distance = list(type = "least-cost", 
                				resistance = resistanceFocal, 
                				mask =  habitatPatchesFocal),
                	 attribute = NULL,
                	 metric = "PC", 
                	 probability = prob, 
                	 distance_thresholds = maxdist,
                	 rasterparallel = T)
     
  # Also transform and prep flux and connector layers for final report figures
flux <- PCfocal[[4]] %>%
		calc(., fun=log)
flux[!is.finite(flux)] <- NA
flux <-	calc(flux, fun=rescaleR)    
 
connector <- PCfocal[[5]] %>%
		calc(., fun=log)
connector[!is.finite(connector)] <- NA
connector <- calc(connector, fun=rescaleR)    
     
                	 
## Save output files
  # Focal area
  #geotif
writeRaster(PCfocal, 
				file.path(outDir, 
				paste0(species, "_PC_FocalArea.tif")), 
				overwrite=TRUE)  

writeRaster(flux, 
				file.path(outDir, 
				paste0(species, "_PC_FocalArea_Flux_log&0-1.tif")), 
				overwrite=TRUE)  

writeRaster(connector, 
				file.path(outDir, 
				paste0(species, "_PC_FocalArea_Connector_log&0-1.tif")), 
				overwrite=TRUE)  		

rm(PCfocal)

} #end loop

# End script


