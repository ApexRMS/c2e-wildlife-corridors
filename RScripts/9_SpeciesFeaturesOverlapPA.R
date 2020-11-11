#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping       
# Calculate overlap between species measures and protected areas 
# 11-2020                                       					
#                             
#	1.Inputs (for focal species):                                   
#    -Patch importance (flux and connector)
#	- habitat suitability
#	-ecopark PA locations
#  
#   Outputs:
#    -Table of overlap values (%)
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################


## Workspace ---------------------------------------------------------

  # Packages
library(tidyverse)
library(raster)
library(sf)

setwd("~/Dropbox/Documents/ApexRMS/Work/A233 - Cootes to Escarpment")

  # Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
outDir <- "Results"

source(file.path("Data/Parameters", "a233_InputParams.R")) # project level parameters

  # EcoParks data
ecopark <- st_read(file.path(
			paste0(rawDataDir, "/Land use land cover/EcoParkLands"), 			"CurrentEcoParkLands.shp")) 
  # Focal area
focalArea <- raster(file.path(procDataDir, "LULC_FocalArea.tif")) %>%
			 calc(., fun=function(x){ifelse(x > 0, 1, x)})


## Load species feature files
  # Habitat suitability, values over 60 only
BLBR_HS <- raster(file.path(procDataDir, 
					  paste0("BLBR", "_HabitatSuitability_FocalArea.tif"))) %>%
			 			calc(., fun=function(x){ifelse(x >= suitabilityThreshold, x, NA)})

EMBL_HS <- raster(file.path(procDataDir, 
					  paste0("EMBL", "_HabitatSuitability_FocalArea.tif"))) %>%
			 			calc(., fun=function(x){ifelse(x >= suitabilityThreshold, x, NA)})
ODVI_HS <- raster(file.path(procDataDir, 
					  paste0("ODVI", "_HabitatSuitability_FocalArea.tif"))) %>%
			 			calc(., fun=function(x){ifelse(x >= suitabilityThreshold, x, NA)})

  # Intra fraction
BLBR_PCflux <- raster(file.path(outDir, 
					  paste0("BLBR", "_PC_FocalArea.tif")),	band = 4) %>%	
					  extend(., extent(focalArea), value=NA) %>%			
					  crop(., focalArea) %>%
					  mask(., focalArea)
EMBL_PCflux <- raster(file.path(outDir, 
					  paste0("EMBL", "_PC_FocalArea.tif")), band = 4) %>%	
					  extend(., extent(focalArea), value=NA) %>%			
					  crop(., focalArea) %>%
					  mask(., focalArea)
ODVI_PCflux <- raster(file.path(outDir, 
					  paste0("ODVI", "_PC_FocalArea.tif")), band = 4) %>%	
					  extend(., extent(focalArea), value=NA) %>%			
					  crop(., focalArea) %>%
					  mask(., focalArea)

  # Connector fraction					  
BLBR_PCconnect <- raster(file.path(outDir, 
					  paste0("BLBR", "_PC_FocalArea.tif")),
					  band = 5) %>%	
					  extend(., extent(focalArea), value=NA) %>%			
					  crop(., focalArea) %>%
					  mask(., focalArea)
EMBL_PCconnect <- raster(file.path(outDir, 
					  paste0("EMBL", "_PC_FocalArea.tif")),
					  band = 5) %>%	
					  extend(., extent(focalArea), value=NA) %>%			
					  crop(., focalArea) %>%
					  mask(., focalArea)
ODVI_PCconnect <- raster(file.path(outDir, 
					  paste0("ODVI", "_PC_FocalArea.tif")),
					  band = 5) %>%	
					  extend(., extent(focalArea), value=NA) %>%			
					  crop(., focalArea) %>%
					  mask(., focalArea)

## Combine parks layer and focal area --------------------------------------------
  # Overlay ecoparks layer
ecoparkSingle <- ecopark %>%
				st_cast(., "POLYGON") %>%
				st_transform(., crs=crs(focalArea))  # Transform to focal CRS
  # Convert to comparable raster
ecoparkRast <- ecoparkSingle %>%
				rasterize(., focalArea, field="OBJECTID") %>%
				crop(., focalArea) %>%
				mask(., focalArea)

  # Create binary version of parks raster
ecoparkBin <- calc(ecoparkRast, fun=function(x){ifelse(x >= 0, 1, x)})
  plot(ecoparkBin)
#  88455/225847 #percentage of focal area

## Calculate  representation ------------------------------------------------------

  # Matrix in which to save representation values, col=species, rows=features
allModelRep <- matrix(NA, nrow=3, ncol=3)
rownames(allModelRep) <- specieslist
colnames(allModelRep) <- c("HabitatSuit", "PIflux", "PIconnect")

  # Per feature/species combination how much occurs in protected areas? 
#BLBR
    #Suitability
  coverage <- overlay(x= ecoparkBin, y = BLBR_HS, fun=function(x, y){(x * y)})
  allModelRep[1, 1] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(BLBR_HS, "sum", na.rm=TRUE)
    # Flux
  coverage <- overlay(x= ecoparkBin, y= BLBR_PCflux, fun=function(x, y){(x * y)})
  allModelRep[1, 2]  <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(BLBR_PCflux, "sum", na.rm=TRUE)
    # Connector
  coverage <- overlay(x= ecoparkBin, y= BLBR_PCconnect, fun=function(x, y){(x * y)})
  allModelRep[1, 3] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(BLBR_PCconnect, "sum", na.rm=TRUE)
#ODVI
    #Suitability
  coverage <- overlay(x= ecoparkBin, y = ODVI_HS, fun=function(x, y){(x * y)})
  allModelRep[2, 1] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(ODVI_HS, "sum", na.rm=TRUE)
    # Flux
  coverage <- overlay(x= ecoparkBin, y= ODVI_PCflux, fun=function(x, y){(x * y)})
  allModelRep[2, 2]  <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(ODVI_PCflux, "sum", na.rm=TRUE)
    # Connector
  coverage <- overlay(x= ecoparkBin, y= ODVI_PCconnect, fun=function(x, y){(x * y)})
  allModelRep[2, 3] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(ODVI_PCconnect, "sum", na.rm=TRUE)
#EMBL
    #Suitability
  coverage <- overlay(x= ecoparkBin, y = EMBL_HS, fun=function(x, y){(x * y)})
  allModelRep[3, 1] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(EMBL_HS, "sum", na.rm=TRUE)
    # Flux
  coverage <- overlay(x= ecoparkBin, y= EMBL_PCflux, fun=function(x, y){(x * y)})
  allModelRep[3, 2]  <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(EMBL_PCflux, "sum", na.rm=TRUE)
    # Connector
  coverage <- overlay(x= ecoparkBin, y= EMBL_PCconnect, fun=function(x, y){(x * y)})
  allModelRep[3, 3] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(EMBL_PCconnect, "sum", na.rm=TRUE)

barplot(t(100*allModelRep), cex.names=0.75, angle=30, beside=TRUE, ylab="% of Partner Lands", xlab="Species", ylim=c(0, 100))
legend("topleft", 
       legend=c("Suitable habitat", "Connectivity (flux)", "Connectivity (stepping stone)"), 
       fill=c("grey10", "grey50", "grey80"))
       
       
## Calculate percentage of focal area that is suitable habitat for each species       
length(which(values(BLBR_HS)>0)) / 225847 #BLBR

length(which(values(ODVI_HS)>0)) / 225847 #ODVI

length(which(values(EMBL_HS)>0)) / 225847 #EMBL

  # End script
