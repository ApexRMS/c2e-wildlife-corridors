#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping            
# Process bridge and culvert data and include in resistance layer                                
# 09/2020
#    
# 1.Process culvert and bridge point data	
#	- create polygon from points 
#	-clip to match focal study area (polygon)
# 2. Calculate resistance values, output crosswalk tables for bridges + culverts
# 3. Convert to raster, buffer
#      
# Inputs:                                                           
#    - C2e Bridges and C2e Culvert attributes and shape files
#	 - Generic species resistance raster 
#	 - area shape files	                                                                                                            
#
# Outputs:     
#	 - OR values                                                     
#    - Bridge and Culvert Resistance layer                       
#                                                                   
# Script created by B Rayfield, C Debeyser, C Tucker for ApexRMS	
#####################################################################

## Workspace ---------------------------------------------------------

  # Packages
library(tidyverse)
library(raster)
library(sf)
library(stringr)

## Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A233 - Cootes to Escarpment" #CT
#projectDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A233_RBGConnectivity/a233" #BR
dataDir <- paste0(projectDir, "/Data/Raw")
outDir <- paste0(projectDir, "/Data/Processed")

## Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

## Input parameters
source(file.path(dataDir, "a233_InputParameters.R")) # project level parameters
	polygonBufferWidth
	roadbuffer
	minCulvertSize
	culvertResistancevalue
	bridgeResistancevalue

## Functions (output units = meters)
  # Attribute: OR calculation for circular culverts. Converts mm measures to m.
ORcirc <- function(length, width){
			(pi * (0.5 * width / 1000)^2) / (length)}
  # Attribute: OR calculation for box culverts
ORbox <- function(length, width, height){
			((height / 1000) * (width / 1000)) / (length)}

## Load data  -------------------------------------------------------

  # Read in generic resistance raster
genericResistanceBuffer <- raster(
							file.path(outDir, "GenericResistance_20km_buffer.tif"))

  # Study area - the polygon that encompases the Cootes to Escarpment EcoPark System
  # Projected file is our desired projection (see 1_StudyAreaLULCL)
focalAreaPolygon <- st_read(
							file.path(outDir, "FocalArea.shp")) # Study area polygon
studyAreaPolygon <- st_read(
							file.path(outDir, "StudyArea_20km.shp")) # Study area polygon

  # Culverts
culvertPoints <- st_read(
						file.path(
						paste0(dataDir, "/Land use land cover/C2E_Culverts"),
						"C2E_Culverts.shp"))
  
  # Bridges
bridgePoints <- st_read(
						file.path(
						paste0(dataDir, "/Land use land cover/C2E_Bridges/"),
                        "C2E_Bridges.shp"))
  
## Adjust bridge and culvert extent and projection to match resistance layer------------

  # Extract bridge length from "Description" column, unit=meters 
bridgePoints$LENGTH <-  bridgePoints$DESCRIPTIO %>% 
  						word(., )  %>%
  						str_remove(., "m") %>%
  						as.numeric(.)

  # Reproject culvert and bridge geo data to same crs as studyAreaPolygon
culvertProjected <- st_transform(culvertPoints, crs(focalAreaPolygon))
bridgeProjected <- st_transform(bridgePoints, crs(focalAreaPolygon))

  # Filter bridge/culvert points to those in focal area only
culvertFocalArea <-   suppressWarnings(st_intersection(culvertProjected, focalAreaPolygon)) # Clip to focal study area
bridgeFocalArea <-   suppressWarnings(st_intersection(bridgeProjected, focalAreaPolygon)) # Clip to focal study area

## Calculate resistance values for culverts -------------------------------------------------------
  
  # Culverts types 
  # "Circular" "Arch"     NA         "Ellipse"  "Box" 
  # Selects between OR_circ and OR_box fct based on type
  # Uses the height and diameter for each culvert 
  	#When 2 measures of height or diameter are available uses the minimum
 
culvertFocalArea$OR <- culvertFocalArea %>%
						apply(., 1, 
							FUN <- function(x){
								ifelse(x$CULVERTSHA!="Box", 
									ORcirc(length = (x$BARREL_LEN), 
										width = (min(x$DIAMETER_W, 
												 x$DIAMETER_1, 
												 na.rm = TRUE))),
									ORbox(length = (x$BARREL_LEN), 
										width = (min(x$DIAMETER_W, x$DIAMETER_1, na.rm=TRUE)), 
										height = min(x$HEIGHT_S1, x$HEIGHT_S2, na.rm=TRUE)))
								})
  
  # About half of the OR values fall less than 0.015. According to Halton Road Ecology paper, values  must be > 0.05 for even smallest mammals    
 
  # Classify all values > minCulvertSize as having resistance = culvertResistancevalue, else remove culvert 
culvertFocalAreaRed <- culvertFocalArea %>% 
						filter(., OR > minCulvertSize) %>%
						mutate(Resistance = culvertResistancevalue)

## Generate crosswalk table for culverts
culvertCross <- culvertFocalAreaRed %>%
					as.data.frame() %>%
					dplyr::select(CULV_ID, ROAD, OR, Resistance) %>%
					filter_all(all_vars(!is.infinite(.)))
				
## Calculate resistance values for bridges ------------------------------------------------------
  # All bridges given resistance value of bridgeResistancevalue
bridgeFocalAreaRed <- bridgeFocalArea %>% 
						mutate(Resistance = bridgeResistancevalue)

 # Generate crosswalk table for bridges
bridgeCross <- bridgeFocalAreaRed %>%
					as.data.frame() %>%
					dplyr::select(TRF_ID, STREET, Resistance) %>%
					filter_all(all_vars(!is.infinite(.)))

## Convert to resistance raster files ------------------------------------------------------

## Create a buffer around culvert and bridge objects.
  # add distance = roadbuffer distance in all directions, to match road buffering
culvertBuffer <- culvertFocalAreaRed %>%
					st_buffer(., roadbuffer)

bridgeBuffer <- bridgeFocalAreaRed %>%
					st_buffer(., roadbuffer)

culvertRaster <- culvertBuffer %>% 
					rasterize(., genericResistanceBuffer, field = "Resistance") # Rasterize

bridgeRaster <- bridgeBuffer %>% 
					rasterize(., genericResistanceBuffer, field = "Resistance") # Rasterize

## Combine rasters by selecting min resistance value per cell
resistanceStack <- stack(genericResistanceBuffer, culvertRaster, bridgeRaster)
combinedRaster  <- min(resistanceStack, na.rm=TRUE)

## Crop to study area --------------------------------------------------------
culvertRasterStudyArea <- culvertRaster %>%
  	crop(., extent(studyAreaPolygon), snap="out") %>% # Crop to focal area extent
  	mask(., mask=studyAreaPolygon) %>% # Clip to focal area
  	trim(.) # Trim extra white spaces

bridgeRasterStudyArea  <- bridgeRaster %>%
  	crop(., extent(studyAreaPolygon), snap="out") %>% # Crop to focal area extent
  	mask(., mask=studyAreaPolygon) %>% # Clip to focal area
  	trim(.) # Trim extra white spaces

combinedRasterStudyArea  <- combinedRaster %>%
  	crop(., extent(studyAreaPolygon), snap="out") %>% # Crop to focal area extent
  	mask(., mask=studyAreaPolygon) %>% # Clip to focal area
  	trim(.) # Trim extra white spaces

## Crop to focal area ----------------------------------------------------------
culvertRasterFocal <- culvertRaster %>%
  	crop(., extent(focalAreaPolygon), snap="out") %>% # Crop to focal area extent
  	mask(., mask=focalAreaPolygon) %>% # Clip to focal area
  	trim(.) # Trim extra white spaces

bridgeRasterFocal  <- bridgeRaster %>%
  	crop(., extent(focalAreaPolygon), snap="out") %>% # Crop to focal area extent
  	mask(., mask=focalAreaPolygon) %>% # Clip to focal area
  	trim(.) # Trim extra white spaces

combinedRasterFocal  <- combinedRaster %>%
  crop(., extent(focalAreaPolygon), snap="out") %>% # Crop to focal area extent
  mask(., mask=focalAreaPolygon) %>% # Clip to focal area
  trim(.) # Trim extra white spaces

## Save outputs ------------------------------------------------------

write.csv(culvertCross, 
				file.path(
				paste0(dataDir, "/Resistance"), 
                "GenericCulvertResistanceCrosswalk.csv"), 
                row.names=F)
write.csv(bridgeCross, 
				file.path(
				paste0(dataDir, "/Resistance"), 
                "GenericBridgeResistanceCrosswalk.csv"), 
                row.names=F)

  # Polygons
st_write(bridgeFocalAreaRed, 
			file.path(outDir, "allBridgesFocal.shp"), 
			driver="ESRI Shapefile")
st_write(culvertFocalArea, 
			file.path(outDir, "allCulvertsFocal.shp"), 
			driver="ESRI Shapefile")  
st_write(culvertFocalAreaRed, 
			file.path(outDir, "suitableCulvertsFocal.shp"), 
			driver="ESRI Shapefile")

  # Rasters of resistance (Study Area buffered) 
writeRaster(culvertRaster, 
			file.path(outDir, "Generic_ResistanceCulvert_20km_buffer.tif"), 
			overwrite=T)
writeRaster(bridgeRaster, 
			file.path(outDir, "Generic_ResistanceBridge_20km_buffer.tif"), 
			overwrite=T)
writeRaster(combinedRaster, 
			file.path(outDir, "Generic_ResistanceCulvertBridge_20km_buffer.tif"), 
			overwrite=T)

  # Rasters of resistance (Study Area)
writeRaster(culvertRasterStudyArea, 
			file.path(outDir, "Generic_ResistanceCulvert_20km.tif"), 
			overwrite=T)
writeRaster(bridgeRasterStudyArea, 
			file.path(outDir, "Generic_ResistanceBridge_20km.tif"), 
			overwrite=T)
writeRaster(combinedRasterStudyArea, 
			file.path(outDir, "Generic_ResistanceCulvertBridge_20km.tif"), 
			overwrite=T)

  # Rasters of resistance (Focal Area)
writeRaster(culvertRasterFocal, 
			file.path(outDir, "Generic_ResistanceCulvert_FocalArea.tif"), 
			overwrite=T)
writeRaster(bridgeRasterFocal, 
			file.path(outDir, "Generic_ResistanceBridge_FocalArea.tif"), 
			overwrite=T)
writeRaster(combinedRasterFocal, 
			file.path(outDir, "Generic_ResistanceCulvertBridge_FocalArea.tif"), 
			overwrite=T)

  # ASCII
    # Rasters of resistance (Study Area buffered)
writeRaster(combinedRaster, 
			file.path(outDir, "Generic_ResistanceCulvertBridge_20km_buffer.asc"), 
			overwrite=T)
