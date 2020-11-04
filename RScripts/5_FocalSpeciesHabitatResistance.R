#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping            
# Create habitat and resistance maps for focal species 
# 09/2020
#                                                                   
# Inputs:                                                           
#    - The combined LULC layer produced by 1_StudyAreaLULC.R        
#    - The LULC/habitat suitability crosswalk          
#	 - Shapefiles for study areas including w buffers 
#	 - Min patch size csv            
# Outputs (for all scales):                                                          
#    - Habitat suitability map                                      
#    - Habitat patch map                                            
#    - Resistance map                                               
#                                                                   
# Script created by Bronwyn Rayfield & C Tucker for ApexRMS                    
#####################################################################

## Workspace ---------------------------------------------------------

  # Packages
library(tidyverse)
library(raster)
library(sf)

  # Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

  # Directories
dataDir <- "Data/Raw"
outDir <- "Data/Processed"

  # Input parameters
source(file.path("Data/Parameters", "a233_InputParams.R")) # project level parameters
	specieslist
	polygonBufferWidth
	suitabilityThreshold

## Load data
  # Combined LULC layers
LULC <- raster(
			file.path(outDir, 
			paste0("LULC_", polygonBufferWidth, "km.tif")))
LULC_buffer <- raster(
				file.path(outDir, 
				paste0("LULC_", polygonBufferWidth, "km_buffered.tif")))

  # Focal area polygon
focalArea <- st_read(file.path(outDir, "FocalArea.shp"))
  # Study area polygon
studyArea <- st_read(file.path(outDir, "studyarea_20km.shp"))

  # Tabular data
crosswalkHabSuit <- read_csv(
					file.path(
					paste0(dataDir, "/Focal Species"), 
					"FocalSpeciesHabitatSuitabilityCrosswalk.csv"))
crosswalkResist <- read_csv(
					file.path(
					paste0(dataDir, "/Resistance"), 
					"FocalSpeciesResistanceCrosswalk.csv"))
minPatchSize <- read_csv(
					file.path(
					paste0(dataDir, "/Focal Species"), 
					"FocalSpeciesMinPatchSize.csv"))

  # Parks
currentParks <- st_read(
				file.path(paste0(dataDir, "/Land use land cover/EcoParkLands"), 
				"CurrentEcoParkLands.shp"))
hendryMulti <- currentParks[which(currentParks$OBJECTID==24), ]
hendryMulti$DestinationID <- 100
hendry <- hendryMulti %>%
			st_cast(., "POLYGON") %>%
		    st_transform(., crs=st_crs(studyArea)) %>% # Transform to SOLRIS CRS
  		    st_intersection(., studyArea) %>%
			rasterize(., LULC_buffer, field = "DestinationID") %>% # Rasterize
			calc(., fun=function(x){ifelse(x==100, 100, NA)})

## Create habitat suitability layer  ---------------------------------------------------------

## For loop over species list to generate focal species habitat suitability and resistance rasters

for(i in specieslist){

species <- i

suitabilityRaster <- LULC_buffer %>%
  					 reclassify(., rcl=crosswalkHabSuit[, c("LULC_ID", species)])

if(species=="EMBL"){
	suitabilityRaster <- max(suitabilityRaster, hendry, na.rm=TRUE)
}


## Create habitat patch layer by species min patch size ---------------------------------

patchSizeThreshold <- minPatchSize$MinPatchSizeHa[minPatchSize$Species == species]

## Generate habitat patch raster
habitatRaster <- Which(suitabilityRaster >= suitabilityThreshold) #habitats above suitability threshold
  #Convert from hectares to m
conversionFromHa <- res(habitatRaster)[1] * res(habitatRaster)[2] * (1 / 10000)
  #Combine neighbouring patches of like suitability
habitatClump <- clump(habitatRaster)
habitatClumpID <- data.frame(freq(habitatClump))
  # Remove clump observations with frequency smaller than minimum habitat patch size (ha)
habitatClumpID <- habitatClumpID[habitatClumpID$count < patchSizeThreshold/conversionFromHa,]
habitatRaster[Which(habitatClump %in% habitatClumpID$value)] <- 0

  # Create raster with ID's for all patches > min threshold & > suitability threshold
habitatRasterCont <- clump(habitatRaster) 

  # Binary raster, suitable habitat
binaryHabitatRaster <- calc(habitatRaster, fun=function(x){ifelse(x > 0, 1, x)})

  
## Create resistance layer ---------------------------------------------------------

  # Reclassify
resistanceRasterReclass <- LULC_buffer %>%
  							reclassify(., rcl=crosswalkResist[, c("LULC_ID", species)])

  # Overlay habitat patches
resistanceRaster <- overlay(resistanceRasterReclass, habitatRaster, 
							fun = function(x, y){return(ifelse(y==1, 1, x))})
#Reclass to assign habitat patches a resistance value = 1 (note that both overlaid values of both 3 and 5 correspond to habitat patches)
#resistanceRaster <- reclassify(resistanceRasterOverlay, rcl = matrix(c(3, 1, 5, 1, 7, 1, 9, 1), ncol=2, byrow = T))

## Crop layers to focal region and study area--------------------------------------------
  
  # Focal Area
suitabilityRasterFocal <- suitabilityRaster %>%
  							crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  							mask(., mask=focalArea) %>% # Clip to focal area
  							trim(.) # Trim extra white spaces  

habitatRasterFocal <- habitatRaster %>%
  						crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
 						mask(., mask=focalArea) %>% # Clip to focal area
 						trim(.) # Trim extra white spaces

habitatRasterContFocal <- habitatRasterCont %>%
  							crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  							mask(., mask=focalArea) %>% # Clip to focal area
  							trim(.) # Trim extra white spaces

binaryHabitatRasterFocal <- binaryHabitatRaster %>%
  							crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  							mask(., mask=focalArea) %>% # Clip to focal area
  							trim(.) # Trim extra white spaces

resistanceRasterFocal <- resistanceRaster %>%
  							crop(., extent(focalArea), snap="out") %>% # Crop to focal area extent
  							mask(., mask=focalArea) %>% # Clip to focal area
  							trim(.) # Trim extra white spaces

  # Study Area
suitabilityRasterStudyArea <- suitabilityRaster %>%
  								crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
  								mask(., mask=studyArea) %>% # Clip to focal area
  								trim(.) # Trim extra white spaces  

habitatRasterStudyArea <- habitatRaster %>%
  							crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
  							mask(., mask=studyArea) %>% # Clip to focal area
  							trim(.) # Trim extra white spaces

habitatRasterContStudyArea <- habitatRasterCont %>%
  								crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
  								mask(., mask=studyArea) %>% # Clip to focal area
  								trim(.) # Trim extra white spaces

resistanceRasterStudyArea <- resistanceRaster %>%
  								crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
  								mask(., mask=studyArea) %>% # Clip to focal area
  								trim(.) # Trim extra white spaces


## Save outputs ---------------------------------------------------------
  # geotif
  # Focal Area
writeRaster(suitabilityRasterFocal, 
			file.path(outDir, 
			paste0(species, "_HabitatSuitability_FocalArea.tif")), 
			overwrite=TRUE)
writeRaster(habitatRasterFocal, 
			file.path(outDir, 
			paste0(species, "_HabitatPatch_FocalArea.tif")), 
			overwrite=TRUE)
writeRaster(binaryHabitatRasterFocal, 
			file.path(outDir, 
			paste0(species, "_HabitatSuitability_Binary_FocalArea.tif")), 
			overwrite=TRUE)			
writeRaster(habitatRasterContFocal, 
			file.path(outDir, 
			paste0(species, "_HabitatPatchID_FocalArea.tif")), 
			overwrite=TRUE)
writeRaster(resistanceRasterFocal, 
			file.path(outDir, 
			paste0(species, "_Resistance_FocalArea.tif")), 
			overwrite=TRUE)

  # Study Area Unbuffered
writeRaster(suitabilityRasterStudyArea, 
			file.path(outDir, 
			paste0(species, "_HabitatSuitability_", polygonBufferWidth, "km.tif")), 
			overwrite=TRUE)
writeRaster(habitatRasterStudyArea, 
			file.path(outDir, 
			paste0(species, "_HabitatPatch_", polygonBufferWidth, "km.tif")), 
			overwrite=TRUE)			
writeRaster(habitatRasterContStudyArea, 
			file.path(outDir, 
			paste0(species, "_HabitatPatchID_",  polygonBufferWidth, "km.tif")), 
			overwrite=TRUE)
writeRaster(resistanceRasterStudyArea, 
			file.path(outDir, 
			paste0(species, "_Resistance_", polygonBufferWidth, "km.tif")), 
			overwrite=TRUE)

  # Study Area Buffered
writeRaster(habitatRaster, 
			file.path(outDir, 
			paste0(species, "_HabitatPatch_", polygonBufferWidth, "km_buffered.tif")), 
			overwrite=TRUE)
writeRaster(resistanceRaster, 
			file.path(outDir, 
			paste0(species, "_Resistance_", polygonBufferWidth, "km_buffered.tif")), 
			overwrite=TRUE)
writeRaster(habitatRasterCont, 
			file.path(outDir, 
			paste0(species, "_HabitatPatchID_", polygonBufferWidth, "km_buffered.tif")), 
			overwrite=TRUE)

} #end loop


## End script