#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping       
# Calculate overall summed important of patches for 3 focal species
# 10-2020                                       					
#                             
#	1.Inputs (for focal species):                                   
#    -Patch importance outputs
#  
#   Outputs:
#    -Raster with Normalized sum score for each patch. 
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

## Workspace ---------------------------------------------------------

  # Packages
library(tidyverse)
library(raster)
library(sf)

  # Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A233 - Cootes to Escarpment" #CT
rawDataDir <- file.path(projectDir, "Data/Raw")
procDataDir <- file.path(projectDir, "Data/Processed")
outDir <- file.path(projectDir, "Results")

  # Functions
rescaleR <- function(x, new.min = 0, new.max = 1) {
   x.min = min(x, na.rm=TRUE)
   x.max = max(x, na.rm=TRUE)
   new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

  # Input parameters
source(file.path(rawDataDir, "a233_InputParameters.R")) # project level parameters
	specieslist

## Load species patch importance files -----------------------------------------

focalArea <- raster(file.path(procDataDir, "LULC_FocalArea.tif"))


  # Habitat suitability
BLBR_HS <- raster(file.path(procDataDir, 
					  paste0("BLBR", "_HabitatSuitability_FocalArea.tif")))
EMBL_HS <- raster(file.path(procDataDir, 
					  paste0("EMBL", "_HabitatSuitability_FocalArea.tif")))
ODVI_HS <- raster(file.path(procDataDir, 
					  paste0("ODVI", "_HabitatSuitability_FocalArea.tif")))

density <- raster(file.path(outDir, "Generic_CulvertBridge_5pairwise_cum_curmap.tif"))

  # Intra fraction
BLBR_PCflux <- raster(file.path(outDir, 
					  paste0("BLBR", "_PC_FocalArea.tif")),
					  band = 4)
EMBL_PCflux <- raster(file.path(outDir, 
					  paste0("EMBL", "_PC_FocalArea.tif")),
					  band = 4)
ODVI_PCflux <- raster(file.path(outDir, 
					  paste0("ODVI", "_PC_FocalArea.tif")),
					  band = 4)

  # Connector fraction					  
BLBR_PCconnect <- raster(file.path(outDir, 
					  paste0("BLBR", "_PC_FocalArea.tif")),
					  band = 5)
EMBL_PCconnect <- raster(file.path(outDir, 
					  paste0("EMBL", "_PC_FocalArea.tif")),
					  band = 5)
ODVI_PCconnect <- raster(file.path(outDir, 
					  paste0("ODVI", "_PC_FocalArea.tif")),
					  band = 5)
					  
## Process raster layers -----------------------------------------

  # Restrict habitat suitability to 60% only & Scale
#
BLBR_HSred <- BLBR_HS %>%
			 calc(., fun=function(x){
			 	ifelse(x >= suitabilityThreshold, x, NA)})
BLBR_HSbin <- calc(BLBR_HSred, fun=function(x){ifelse(x > 0, 1, x)})

BLBR_HSsc <- BLBR_HSred  %>%
			 scale(., center=TRUE, scale=TRUE)	%>%
			 calc(., fun = rescaleR)
#			 
EMBL_HSred <- EMBL_HS  %>%
			 calc(., fun=function(x){
			 	ifelse(x>= suitabilityThreshold, x, NA)}) 
EMBL_HSbin <- calc(EMBL_HSred, fun=function(x){ifelse(x > 0, 1, x)})

EMBL_HSsc <- EMBL_HSred  %>%
			 scale(., center=TRUE, scale=TRUE) %>%
			 calc(., fun = rescaleR)	
#			  
ODVI_HSred <- ODVI_HS %>% 
			 calc(., fun=function(x){
			 	ifelse(x>= suitabilityThreshold, x, NA)}) 
ODVI_HSbin <- calc(ODVI_HSred, fun=function(x){ifelse(x > 0, 1, x)})

ODVI_HSsc <- ODVI_HSred  %>%
			 scale(., center=TRUE, scale=TRUE)	%>%
			 calc(., fun = rescaleR)

## Make combined habitat suitability layer
combinedHS <- sum(stack(BLBR_HSbin, EMBL_HSbin, ODVI_HSbin), na.rm=TRUE)
combinedHS <- calc(combinedHS, fun=function(x){ifelse(x==0, NA, 1)})

## Crop, scale and log density
densityCrop <- density %>%
				crop(., focalArea) %>%
				mask(., combinedHS) %>%  
				calc(., fun=function(x){log(x)}) %>% #log values for normality
				scale(., center=TRUE, scale=TRUE) %>%
			 	calc(., fun = rescaleR)	

##Crop, scale and log patch importance values
BLBR_PCfluxExt <-  BLBR_PCflux %>%
					extend(., extent(focalArea), value=NA) %>%
					scale(., center=TRUE, scale=TRUE) %>%
			 		calc(., fun = rescaleR)
					
EMBL_PCfluxExt <-  EMBL_PCflux %>%
					extend(., extent(focalArea), value=NA) %>%
					scale(., center=TRUE, scale=TRUE) %>%
			 		calc(., fun = rescaleR)
					
ODVI_PCfluxExt <-  ODVI_PCflux %>%
					extend(., extent(focalArea), value=NA) %>%
					scale(., center=TRUE, scale=TRUE)%>%
			 		calc(., fun = rescaleR)
					
  #					
BLBR_PCconnectExt <-  BLBR_PCconnect %>%
						extend(., extent(focalArea), value=NA) %>%
						scale(., center=TRUE, scale=TRUE) %>%
			 			calc(., fun = rescaleR)
					
EMBL_PCconnectExt  <-  EMBL_PCconnect %>%
						extend(., extent(focalArea), value=NA) %>%
						scale(., center=TRUE, scale=TRUE) %>%
			 			calc(., fun = rescaleR)
					
ODVI_PCconnectExt  <-  ODVI_PCconnect %>%
						extend(., extent(focalArea), value=NA) %>%
						scale(., center=TRUE, scale=TRUE) %>%
			 			calc(., fun = rescaleR)
	
	# Ignore warnings

## Combine normalized data layers and calculate the sum of all layers ----------------

  # Combine all species
allSp <- stack(BLBR_PCfluxExt, BLBR_PCconnectExt, EMBL_PCfluxExt, EMBL_PCconnectExt, ODVI_PCfluxExt, ODVI_PCconnectExt, BLBR_HSsc, EMBL_HSsc, ODVI_HSsc, densityCrop)
plot(allSp)
  # Habitat suitability
HSlayers <- stack(BLBR_HSsc, EMBL_HSsc, ODVI_HSsc)
  # Patch importance layer 
PIlayers <- stack(BLBR_PCfluxExt, BLBR_PCconnectExt, EMBL_PCfluxExt, EMBL_PCconnectExt, ODVI_PCfluxExt, ODVI_PCconnectExt)

## Calculate sum layers for each

  # All species 
allSp_PCRawSum <- sum(allSp, na.rm=TRUE)
allSp_SumCrop <- allSp_PCRawSum %>%
				crop(., focalArea) %>%  
				mask(., combinedHS)
  # Range scale - 10 layers, each scaled from 0-1
allSp_range <- calc((allSp_SumCrop/10), fun = rescaleR)
plot(allSp_range)

  #Habitat suitability
HS_PCRawSum <- sum(HSlayers, na.rm=TRUE)
HS_SumCrop <- HS_PCRawSum %>%
				crop(., focalArea) %>%
				mask(., combinedHS)
  # Range scale - 10 layers, each scaled from 0-1
HS_range <- calc((HS_SumCrop/3), fun = rescaleR)
plot(HS_range)

  #Patch importance
PI_PCRawSum <- sum(PIlayers, na.rm=TRUE)
PI_SumCrop <- PI_PCRawSum %>%
				crop(., focalArea) %>%
				mask(., combinedHS)
  # Range scale - 10 layers, each scaled from 0-1
PI_range <- calc((PI_SumCrop/6), fun = rescaleR)
plot(PI_range)


## Save final allSp and allSpRange raster layer-------------------------------------

#writeRaster(allSp_SumCrop, 
#			file.path(outDir, "allSp_SyntheticRaw.tif"), 
#			overwrite=TRUE)
writeRaster(allSp_range, 
			file.path(outDir, "allSp_SyntheticRange.tif"), 
			overwrite=TRUE)			
writeRaster(HS_range, 
			file.path(outDir, "HabitatSuitability_SyntheticRange.tif"), 
			overwrite=TRUE)		
writeRaster(PI_range, 
			file.path(outDir, "PatchImportance_SyntheticRange.tif"), 
			overwrite=TRUE)						