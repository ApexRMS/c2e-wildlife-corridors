#####################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping            
# Process Bridge and Culvert data                                   
# # Script created by B Rayfield, C Debeyser, C Tucker
#    
# 1.Process culvert and bridge point data	
#	- create polygon from points 
#	-clip to match focal study area (polygon)
# 2. Calculate resistance values, output crosswalk tables for bridges + culverts
# 3. Convert to raster, buffer
# 4. Output updated resistance paper
#                         
# Inputs:                                                           
#    - C2e Bridges and C2e Bridges attributes and shape files                                                                                                            #
# Outputs:     
#	 - OR values                                                     
#    - Bridge and Culvert Resistance layer                       
#                                                                   
# for ApexRMS	09/2020
#####################################################################

# Workspace ---------------------------------------------------------

# Packages
library(tidyverse)
library(raster)
library(sf)
library(stringr)

## Set up -----
	# Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A233 - Cootes to Escarpment/"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)


# Functions (output units = meters)
  	# Attribute: OR calculation for circular culverts. Converts mm measures to m.
ORcirc <- function(length, width){(pi * (0.5*width/1000)^2) / (length)}
  	# Attribute: OR calculation for box culverts
ORbox <- function(length, width, height){((height/1000) * (width/1000)) / (length)}


# Read in raster output example to have for resolution details (probably not necessary?)
exRaster <- raster(paste0(outDir, "/Resistance_buffer20km_buffered.tif"))


## Read in data & process -----------
  	# Study area - the polygon that encompases the Cootes to Escarpment EcoPark System
  	#Projected file is our desired projection (see 1_StudyAreaLULCL)
focalAreaPolygon <- st_read(file.path(outDir, "FocalArea.shp")) # Study area polygon
  
  	#Resistance crosswalk
  #crosswalk <- read_csv(file.path(dataDir, "ResistanceCrosswalk.csv"))
  
  	#Culverts
culvertPoints <- st_read(paste0(file.path(dataDir, "C2E_Culverts/"), "C2E_Culverts.shp"))
  
  	#Bridges
bridgePoints <- st_read(paste0(file.path(dataDir, "C2E_Bridges/"), "C2E_Bridges.shp"))
  
  	#Extract bridge length from "Description" column, unit=meters 
bridgePoints$LENGTH <-  bridgePoints$DESCRIPTIO %>% 
  						word(., )  %>%
  						str_remove(., "m") %>%
  						as.numeric(.)
  	#log distribution of length, most <20m

# Reproject culvert and bridge geo data to match focal area  
  	# Reproject to same crs as studyAreaPolygon
culvertProjected <- st_transform(culvertPoints, crs(focalAreaPolygon))
bridgeProjected <- st_transform(bridgePoints, crs(focalAreaPolygon))

# Filter bridge/culvert points to those in focal area only
  	#Select bounding  box 
culvertFocalArea <-   st_intersection(culvertProjected, focalAreaPolygon) # Clip to focal study area
bridgeFocalArea <-   st_intersection(bridgeProjected, focalAreaPolygon) # Clip to focal study area



## Calculate resistance values for culverts ----------
  
  # Culverts types 
  # "Circular" "Arch"     NA         "Ellipse"  "Box" 
  # Default all but "Box" to OR_circ fct
  # Take the minimum height and diameter for each culvert when 2 measures available
 
culvertFocalArea$OR <- culvertFocalArea %>%
						apply(., 1, 
						FUN <- function(x){ifelse(x$CULVERTSHA!="Box", 
						ORcirc(length = (x$BARREL_LEN), 
						width=(min(x$DIAMETER_W, 
						x$DIAMETER_1, 
						na.rm=TRUE))),
						ORbox(length = (x$BARREL_LEN), 
						width=(min(x$DIAMETER_W, x$DIAMETER_1, na.rm=TRUE)), 
						height= min(x$HEIGHT_S1, x$HEIGHT_S2, na.rm=TRUE)))})
  
  # NaNs or Inf for L18, 30, 68, and 97 due to 0 values or NA values in data. Best approach?
  # About half of the OR values fall less than 0.015. According to Halton Road Ecology paper, values  must be > 0.05 for even smallest mammal and Herpetofauna    
  # Classify all values > 0.05 == 10, else remove culvert 

culvertFocalAreaRed <- culvertFocalArea %>% 
					filter(., OR > 0.05) %>%
					mutate(Resistance = 10)


# Crosswalk table for culverts
culvertCross <- as.data.frame(culvertFocalAreaRed)[, c("CULV_ID", "ROAD", "OR", "Resistance")]
write.csv(culvertCross[complete.cases(culvertCross),], file = file.path(outDir, "/CulvertResistanceCrosswalk.csv"), row.names=F)

## Calculate resistance values for bridges ----------
  # All bridges given resistance value of 10
  bridgeFocalAreaRed <- bridgeFocalArea %>% 
						mutate(Resistance = 10)

 # Crosswalk table for bridges
bridgeCross <- as.data.frame(bridgeFocalAreaRed)[, c("TRF_ID", "STREET", "TYPE", "LENGTH", "Resistance")]
write.csv(bridgeCross, file = file.path(outDir, "/BridgeResistanceCrosswalk.csv"), row.names=F)

## Convert to raster files ----------

# Create a buffer around culvert and bridge objects.
  # 12m in all directions, to roughly match road buffering
culvterBuffer <- culvertFocalAreaRed %>%
				st_buffer(., 12)

bridgeBuffer <- bridgeFocalAreaRed %>%
				st_buffer(., 12)


culvertRaster <- culvterBuffer %>% 
				rasterize(., exRaster, field="Resistance") # Rasterize

bridgeRaster <- bridgeBuffer %>% 
				rasterize(., exRaster, field="Resistance") # Rasterize
				

# Save outputs --------------------------
  # Polygons
st_write(bridgeFocalAreaRed, file.path(outDir, "allBridges.shp"), driver="ESRI Shapefile", overwrite=T, append=F)
st_write(culvertFocalArea, file.path(outDir, "allCulverts.shp"), driver="ESRI Shapefile", overwrite=T, append=F)  
st_write(culvertFocalAreaRed, file.path(outDir, "suitableCulverts.shp"), driver="ESRI Shapefile", overwrite=T, append=F)
  # Rasters of resistance
writeRaster(culvertRaster, file.path(outDir, "culvertResistance_20km_buffered.tif"), overwrite=T)
writeRaster(bridgeRaster, file.path(outDir, "bridgeResistance_20km_buffered.tif"), overwrite=T)
