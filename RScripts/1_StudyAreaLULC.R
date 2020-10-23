#####################################################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping                                             
# Create study area from series of points and use it to create one LULC map                          
# 08/2020  
#                                                                                              
# 	1. Create study area                                                                               
# 		- create polygon from points                                                                       
# 		- Buffer polygon by 20km for study area                                                            
# 		- Buffer study area by an additional 20% for connectivity analysis (to remove edge effects)        
# 
# 	2. LULC crop and standardize rasters                                                               
# 		- process all LULC layers: SOLRIS, OLCDB, ORN                                                      
# 		- clip all LULC layers to buffered study area and reclassify                                       
# 	3. Combine all LULC rasters                                                                        
# 	4. Clip all intermediate and final rasters to study area for visualizations                        #                                                                                                  
# Script created by Bronwyn Rayfield, Chloé Debeyser, Caroline Tucker for ApexRMS                                  
#####################################################################################################


## Workspace -------------------------------------

  # Packages
library(tidyverse)
library(sf)
library(raster)


  # Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"

  # Input parameters
source(file.path(rawDataDir, "a233_InputParameters.R"))  #load project level parameters
	polygonBufferWidth
	
## Load data
  # Study area
studyAreaPoints <- read_csv(
				file.path(
				paste0(rawDataDir, "/Study Area"), 
				"Vertices 20200610 DAG.csv")) # Study area vertices
  # LULC data
SOLRISRaw <- raster(
				file.path(
				paste0(rawDataDir,"/Land use land cover/SOLRIS_Version_3_0"), 
				"SOLRIS_Version_3_0_LAMBERT.tif"))
OLCDBRaw <- raster(
				file.path(
				paste0(rawDataDir, "/Land use land cover/Provincial Land Cover 1996 - 28 Class Grid"), 
				"plc1996_28.tif"))
ORNRaw <- st_read(
				file.path(
				paste0(rawDataDir, "/Land use land cover/Ontario_Road_Network__ORN__Segment_With_Address-shp"), 
				"Ontario_Road_Network__ORN__Segment_With_Address.shp"))
OHNRaw <- st_read(
				file.path(
				paste0(rawDataDir, "/Land use land cover/Ontario_Hydro_Network__OHN__-_Waterbody-shp"),
				"Ontario_Hydro_Network__OHN__-_Waterbody.shp"))
urbanRaw <- st_read(
				file.path(
				paste0(rawDataDir, "/Land use land cover/Built-Up_Area-shp"), 
				"Built-Up_Area.shp"))

railRaw <- st_read(
				file.path(
				paste0(rawDataDir, "/Land use land cover/ORWNTRK/LIO-2019-09-30"), 
				"ORWN_TRACK.shp"))

## Additional files
currentParks <- st_read(
				file.path(paste0(rawDataDir, "/Land use land cover/EcoParkLands"), 
				"CurrentEcoParkLands.shp"))
 
  # Hopkins tract + Berry's tract
berryII <- currentParks[which(currentParks$Name_1=="Berry Tract II"),]
hopkinsMulti <- currentParks[which(currentParks$Name_1=="Hopkins Track"),]
hopkins <- st_cast(hopkinsMulti, "POLYGON")
  
  # Resistance crosswalk
crosswalk <- read_csv(
				file.path(
				paste0(rawDataDir, "/Resistance"), 
				"GenericResistanceCrosswalk.csv"))

## Generate focal study area -------------------------------------

  # Create spatial polygon from points
polygon <- studyAreaPoints %>%
  		st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  		summarise(geometry = st_combine(geometry)) %>%
  		st_cast("POLYGON")

  # Reproject polygon to the same crs as lulcRaw
polygonProjected <- st_transform(polygon, crs(SOLRISRaw))

  # Create buffer around polygon for the study area using polygon buffer width
studyArea <- st_buffer(polygonProjected, (polygonBufferWidth * 1000))

  # Create bounding box around study area
	# Calculate 20% of the average of the max length and width of the bounding box to be the study area buffer width
	# Connectivity analyses will be run in this buffered study area to minimize edge effects
studyAreaBBox <- st_bbox(studyArea)
studyAreaBufferWidth <- round(((studyAreaBBox$ymax - studyAreaBBox$ymin) + (studyAreaBBox$xmax - studyAreaBBox$xmin)) / 2 * 0.2, digits = 0) # in m
studyAreaBuffer <- st_buffer(studyArea, studyAreaBufferWidth)

## Crop LULC data to STUDY AREA BUFFER region and standardize rasters ---------------------

  # SOLRIS data
  # Crop to buffered study area
SOLRIS_buffer <- SOLRISRaw %>%
  		crop(., extent(studyAreaBuffer), snap="out") %>% # Crop SOLRIS to study area buffer extent
  		mask(., mask=studyAreaBuffer) %>% # Clip to buffered study area
  		trim(.) # Trim extra white spaces

  # Reclassify values resistance scores
SOLRIS_rcl_buffer <- SOLRIS_buffer %>%
  		reclassify(., rcl=crosswalk[which(crosswalk$Source_Name == "SOLRIS"), c("Source_ID", "Destination_ID")])

  # OLCDB
  # Align with SOLRIS and buffer SOLRIS extent
ext_OLCDB <- extend(extent(SOLRIS_buffer), 1000)

  # Crop to buffered study area
OLCDB_buffer <- OLCDBRaw %>% 
  		crop(., y=ext_OLCDB, snap="out") %>% # Crop OLCDB to new extent 
  		projectRaster(from=., crs=crs(SOLRIS_buffer), method='ngb') %>% # Project to SOLRIS CRS 
  		resample(., SOLRIS_buffer, method='ngb') %>% # Resample to SOLRIS resolution
  		mask(., mask=studyAreaBuffer) %>% # Clip to buffered study area
  		trim(.) # Trim extra white spaces

  # Reclassify
OLCDB_rcl_buffer <- OLCDB_buffer %>%
  		reclassify(., rcl=crosswalk[which(crosswalk$Source_Name == "OLCDB"), c("Source_ID", "Destination_ID")])

  # ORN
   # Crop to buffered study area
ORN_buffer <- ORNRaw %>%
  		st_transform(., crs=st_crs(SOLRIS_buffer)) %>% # Transform to SOLRIS CRS
  		st_intersection(., studyAreaBuffer) %>% # Clip to buffered study area
  		st_buffer(., roadbuffer) # Buffer road lines by constant so that rasterizing aligns with roads in SOLRIS

ORN_bufferMajor <- ORN_buffer[ORN_buffer$ROAD_CLASS %in% c("Arterial", "Collector", "Freeway", "Expressway / Highway", "Ramp", "Rapid Transit"),]

ORN_bufferMinor <- ORN_buffer[ORN_buffer$ROAD_CLASS %in% c("Local / Strata", "Local / Street", "Service", "Alleyway / Laneway", "<NA>", "Local / Unknown"),]

  		
  # Reclassify
ORN_rcl_buffer <- ORN_buffer %>%
  		left_join(., crosswalk[which(crosswalk$Source_Name == "ORN"), c("Destination_ID", "Source_Label")], by=c("ROAD_CLASS" = "Source_Label")) %>% # Add Destination_ID
  		rasterize(., SOLRIS_buffer, field="Destination_ID") # Rasterize

ORN_rcl_bufferMajor <- ORN_bufferMajor %>%
  		left_join(., crosswalk[which(crosswalk$Source_Name == "ORN"), c("Destination_ID", "Source_Label")], by=c("ROAD_CLASS" = "Source_Label")) %>% # Add Destination_ID
  		rasterize(., SOLRIS_buffer, field="Destination_ID") # Rasterize

ORN_rcl_bufferMinor <- ORN_bufferMinor %>%
  		left_join(., crosswalk[which(crosswalk$Source_Name == "ORN"), c("Destination_ID", "Source_Label")], by=c("ROAD_CLASS" = "Source_Label")) %>% # Add Destination_ID
  		rasterize(., SOLRIS_buffer, field="Destination_ID") # Rasterize

  # Railway
   # Crop to buffered study area
rail_buffer <- railRaw %>%
  		st_transform(., crs=st_crs(SOLRIS_buffer)) %>% # Transform to SOLRIS CRS
  		st_intersection(., studyAreaBuffer) # Clip to buffered study area
rail_buffer$TYPE <- "Railway"
rail_rcl_buffer <- rail_buffer %>%
  		left_join(., crosswalk[which(crosswalk$Source_Name == "ORWNTRK"), c("Destination_ID", "Source_Label")], by=c("TYPE" = "Source_Label")) %>% # Add Destination_ID
  		rasterize(., SOLRIS_buffer, field="Destination_ID") # Rasterize  		

  # OHN
   # Crop to buffered study area
OHN_buffer <- OHNRaw %>%
  		st_transform(., crs=st_crs(SOLRIS_buffer)) %>% # Transform to SOLRIS CRS
  		mutate(AREA = as.numeric(st_area(.))) %>% # Calculate the area of each waterbody
  		dplyr::filter(AREA >= 100) %>% # Keep waterbodies larger than 100m2
  		st_intersection(., studyAreaBuffer) # Clip to buffered study area

  # Reclassify
OHN_rcl_buffer <- OHN_buffer %>%
  		left_join(., crosswalk[which(crosswalk$Source_Name == "OHN"), c("Destination_ID", "Source_Label")], by=c("WATERBODY_" = "Source_Label")) %>% # Add Destination_ID
  		rasterize(., SOLRIS_buffer, field="Destination_ID") # Rasterize

  # Urban Areas
urban_buffer <- urbanRaw %>%
  		st_transform(., crs=st_crs(SOLRIS_buffer)) %>% # Transform to SOLRIS CRS
      st_buffer(., dist = 0) %>%
      st_intersection(., studyAreaBuffer)

  # Reclassify
urban_rcl_buffer <- urban_buffer %>%
  		left_join(., crosswalk[which(crosswalk$Source_Name == "Built-up Areas"), c("Destination_ID", "Source_Label")], by=c("COMMUNIT_1" = "Source_Label")) %>% # Add Destination_ID
  		rasterize(., SOLRIS_buffer, field="Destination_ID") # Rasterize

berry_buffer <- berryII %>%
		st_transform(., crs=st_crs(SOLRIS_buffer)) %>% # Transform to SOLRIS CRS
  		st_intersection(., studyAreaBuffer)
berry_buffer$DestinationID <- 194  #Pasture and Abandoned fields		
berry_rcl_buffer <- rasterize(berry_buffer, SOLRIS_buffer, field = "DestinationID") # Rasterize

hopkins_buffer <- hopkins %>%
		st_transform(., crs=st_crs(SOLRIS_buffer)) %>% # Transform to SOLRIS CRS
  		st_intersection(., studyAreaBuffer)
hopkins_buffer$DestinationID <- 194  #Pasture and Abandoned fields		
hopkins_rcl_buffer <- rasterize(hopkins_buffer, SOLRIS_buffer, field = "DestinationID") # Rasterize


## LULC - Combine rasters ---------------------------

#load(file="LULCtemp_line187.RData")
  # Remove undifferentiated cells from SOLRIS so that they can take the value from OLCDB
SOLRIS_rcl_buffer[SOLRIS_rcl_buffer == 250] <- NA

  # Merge rasters: urban gets priority, followed by ORN (Major roads), ORN (Minor roads), OHN, SOLRIS, and last OLCDB
#LULC_buffer <- merge(urban_rcl_buffer, rail_rcl_buffer, ORN_rcl_bufferMajor,  ORN_rcl_bufferMinor, OHN_rcl_buffer, SOLRIS_rcl_buffer, OLCDB_rcl_buffer)
LULC_buffer <- merge(urban_rcl_buffer, rail_rcl_buffer, ORN_rcl_bufferMajor,  ORN_rcl_bufferMinor, OHN_rcl_buffer, berry_rcl_buffer, hopkins_rcl_buffer, SOLRIS_rcl_buffer, OLCDB_rcl_buffer)

## LULC - Order for visualization
LULC_bufferViz <- merge(rail_rcl_buffer, ORN_rcl_bufferMajor, ORN_rcl_bufferMinor, OHN_rcl_buffer, berry_rcl_buffer, hopkins_rcl_buffer, SOLRIS_rcl_buffer, OLCDB_rcl_buffer)

## Crop data to STUDY AREA ------------------------------

  # SOLRIS
SOLRIS <- SOLRISRaw %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

SOLRIS_rcl <- SOLRIS_rcl_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # OLCDB
ext_OLCDB <- extend(extent(SOLRIS), 1000)
  # Crop to buffered study area

OLCDB <- OLCDBRaw %>% 
  crop(., y=ext_OLCDB, snap="out") %>% # Crop OLCDB to new extent 
  projectRaster(from=., crs=crs(SOLRIS), method='ngb') %>% # Project to SOLRIS CRS 
  resample(., SOLRIS, method='ngb') %>% # Resample to SOLRIS resolution
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

OLCDB_rcl <- OLCDB_rcl_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # OHN
OHN <- OHN_buffer %>%
  st_intersection(., studyArea) # Clip to study area

OHN_rcl <- OHN_rcl_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # ORN
ORN <- ORN_buffer %>%
  st_intersection(., studyArea) # Clip to study area

ORN_rcl <- ORN_rcl_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # ORN Major
ORNMajor <- ORN_bufferMajor %>%
  st_intersection(., studyArea) # Clip to study area

ORN_rclMajor <- ORN_rcl_bufferMajor %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # ORN Minor
ORNMinor <- ORN_bufferMinor %>%
  st_intersection(., studyArea) # Clip to study area

ORN_rclMinor <- ORN_rcl_bufferMinor %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # Urban
urban <- urban_buffer %>%
  st_intersection(., studyArea) # Clip to study area

urban_rcl <- urban_rcl_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # Rail
rail <- rail_buffer %>%
  st_intersection(., studyArea) # Clip to study area

rail_rcl <- rail_rcl_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # LULC
LULC <- LULC_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.)

  # LULC viz
LULC_viz <- LULC_bufferViz %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.)


## Crop data to FOCAL AREA ------------------------------

  # SOLRIS
SOLRIS_Focal <- SOLRISRaw %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces

SOLRIS_rcl_Focal <- SOLRIS_rcl_buffer %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # OLCDB
OLCDB_Focal <- OLCDBRaw %>% 
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces

OLCDB_rcl_Focal <- OLCDB_rcl_buffer %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # OHN
OHN_Focal <- OHN_buffer %>%
  st_intersection(., polygonProjected) # Clip to study area

OHN_rcl_Focal <- OHN_rcl_buffer %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # ORN
ORN_Focal <- ORN_buffer %>%
  st_intersection(., polygonProjected) # Clip to study area

ORN_rcl_Focal <- ORN_rcl_buffer %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces
  
  # ORN Major
ORNMajor_Focal <- ORN_bufferMajor %>%
  st_intersection(., polygonProjected) # Clip to study area

ORN_rclMajor_Focal <- ORN_rcl_bufferMajor %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # ORN Minor
ORNMinor_Focal <- ORN_bufferMinor %>%
  st_intersection(., polygonProjected) # Clip to study area

ORN_rclMinor_Focal <- ORN_rcl_bufferMinor %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces
 
  # Urban
urban_Focal <- urban_buffer %>%
  st_intersection(., polygonProjected) # Clip to study area

urban_rcl_Focal <- urban_rcl_buffer %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces

  # Rail
rail_Focal <- rail_buffer %>%
  st_intersection(., polygonProjected) # Clip to study area

rail_rcl_Focal <- rail_rcl_buffer %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.) # Trim extra white spaces
  
  # LULC
LULC_Focal <- LULC_buffer %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.)

LULC_vizFocal <- LULC_bufferViz %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop to study area extent
  mask(., mask= polygonProjected) %>% # Clip to study area
  trim(.)

## Calculate percentage of each class in LULC maps
freqLULC_viz <- as_tibble(freq(LULC_viz)) %>%
  filter(!is.na(value)) %>% mutate(percent = count/sum(count)*100) %>% 
  arrange(desc(percent))

freqLULC_vizFocal <- as_tibble(freq(LULC_vizFocal)) %>%
  filter(!is.na(value)) %>% mutate(percent = count/sum(count)*100) %>% 
  arrange(desc(percent))


## Save outputs --------------------------
  # Tabular
write_csv(freqLULC_viz, 
          file.path(procDataDir, "LULC_Visualization_20km_Freq.csv"))

write_csv(freqLULC_vizFocal, 
          file.path(procDataDir, "LULC_Visualization_FocalArea_Freq.csv"))

  # Save intermediate outputs
st_write(polygonProjected, 
				file.path(procDataDir, "FocalArea.shp"), 
				driver="ESRI Shapefile",
				append = FALSE)
st_write(studyArea, 
				file.path(procDataDir, 
				paste0("StudyArea_", polygonBufferWidth, "km.shp")), 
				driver="ESRI Shapefile",
				append = FALSE)
st_write(studyAreaBuffer, 
				file.path(procDataDir, 
				paste0("StudyArea_", polygonBufferWidth, "km_buffered.shp")), 
				driver="ESRI Shapefile",
				append = FALSE)

  # Focal Area
writeRaster(SOLRIS_Focal, 
				file.path(procDataDir, "SOLRIS_FocalArea.tif"),
				overwrite = TRUE)
writeRaster(SOLRIS_rcl_Focal, 
				file.path(procDataDir,"SOLRIS_reclass_FocalArea.tif"),
				overwrite = TRUE)
writeRaster(OLCDB_Focal, 
				file.path(procDataDir, "OLCDB_FocalArea.tif"),
				overwrite = TRUE)
writeRaster(OLCDB_rcl_Focal, 
				file.path(procDataDir, "OLCDB_reclass_FocalArea.tif"),
				overwrite = TRUE)
writeRaster(ORN_rcl_Focal, 
				file.path(procDataDir, "ORN_reclass_FocalArea.tif"),
				overwrite = TRUE)
writeRaster(ORN_rclMajor, 
				file.path(procDataDir, "ORNMajor_reclass_FocalArea.tif"),
				overwrite = TRUE)	
writeRaster(ORN_rclMinor, 
				file.path(procDataDir, "ORNMinor_reclass_FocalArea.tif"),
				overwrite = TRUE)								
writeRaster(OHN_rcl_Focal, 
				file.path(procDataDir, "OHN_reclass_FocalArea.tif"),
				overwrite = TRUE)
writeRaster(urban_rcl_Focal, 
				file.path(procDataDir, "Urban_reclass_FocalArea.tif"),
				overwrite = TRUE)
writeRaster(rail_rcl_Focal, 
				file.path(procDataDir, "Rail_reclass_FocalArea.tif"),
				overwrite = TRUE)


  # ESRI shapefile				
st_write(ORN_Focal, 
				file.path(procDataDir,"ORN_FocalArea.shp"), 
				driver="ESRI Shapefile",
				append = TRUE)
st_write(ORNMajor_Focal, 
				file.path(procDataDir,"ORNMajor_FocalArea.shp"), 
				driver="ESRI Shapefile",
				append = TRUE)
st_write(ORNMinor_Focal, 
				file.path(procDataDir,"ORNMinor_FocalArea.shp"), 
				driver="ESRI Shapefile",
				append = TRUE)								
st_write(OHN_Focal, 
				file.path(procDataDir, "OHN_FocalArea.shp"), 
				driver="ESRI Shapefile",
				append = TRUE)
st_write(urban_Focal, 
				file.path(procDataDir, "Urban_FocalArea.shp"), 
				driver="ESRI Shapefile",
				append = TRUE)
st_write(rail_Focal, 
				file.path(procDataDir, "Rail_FocalArea.shp"), 
				driver="ESRI Shapefile",
				append = TRUE)


  # Study Area
writeRaster(SOLRIS, 
				file.path(procDataDir, 
				paste0("SOLRIS_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)
writeRaster(SOLRIS_rcl, 
				file.path(procDataDir, 
				paste0("SOLRIS_reclass_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)
writeRaster(OLCDB, 
				file.path(procDataDir, 
				paste0("OLCDB_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)
writeRaster(OLCDB_rcl, 
				file.path(procDataDir, 
				paste0("OLCDB_reclass_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)
writeRaster(ORN_rcl, 
				file.path(procDataDir, 
				paste0("ORN_reclass_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)
writeRaster(ORN_rclMajor, 
				file.path(procDataDir, 
				paste0("ORNMajor_reclass_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)
writeRaster(ORN_rclMinor, 
				file.path(procDataDir, 
				paste0("ORNMinor_reclass_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)				
writeRaster(OHN_rcl, 
				file.path(procDataDir, 
				paste0("OHN_reclass_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)
writeRaster(urban_rcl, 
				file.path(procDataDir, 
				paste0("Urban_reclass_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)
writeRaster(rail_rcl, 
				file.path(procDataDir, 
				paste0("Rail_reclass_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)

  # ESRI shapefile				
st_write(ORN, 
				file.path(procDataDir, 
				paste0("ORN_", polygonBufferWidth, "km.shp")), 
				driver="ESRI Shapefile",
				append = TRUE)
st_write(ORNMajor, 
				file.path(procDataDir, 
				paste0("ORNMajor_", polygonBufferWidth, "km.shp")), 
				driver="ESRI Shapefile",
				append = TRUE)
st_write(ORNMinor, 
				file.path(procDataDir, 
				paste0("ORNMinor_", polygonBufferWidth, "km.shp")), 
				driver="ESRI Shapefile",
				append = TRUE)								
st_write(OHN, 
				file.path(procDataDir, 
				paste0("OHN_", polygonBufferWidth, "km.shp")), 
				driver="ESRI Shapefile",
				append = TRUE)
st_write(urban, 
				file.path(procDataDir, 
				paste0("Urban_", polygonBufferWidth, "km.shp")), 
				driver="ESRI Shapefile",
				append = TRUE)
st_write(rail, 
				file.path(procDataDir, 
				paste0("Rail_", polygonBufferWidth, "km.shp")), 
				driver="ESRI Shapefile",
				append = TRUE)

  # Study Area Buffer
writeRaster(SOLRIS_buffer, 
				file.path(procDataDir, 
				paste0("SOLRIS_", polygonBufferWidth, "km_buffered.tif")),
				overwrite = TRUE)
writeRaster(SOLRIS_rcl_buffer, 
				file.path(procDataDir, 
				paste0("SOLRIS_reclass_", polygonBufferWidth, "km_buffered.tif")),
				overwrite = TRUE)
writeRaster(OLCDB_buffer, 
				file.path(procDataDir, 
				paste0("OLCDB_", polygonBufferWidth, "km_buffered.tif")),
				overwrite = TRUE)
writeRaster(OLCDB_rcl_buffer, 
				file.path(procDataDir, 
				paste0("OLCDB_reclass_", polygonBufferWidth, "km_buffered.tif")),
				overwrite = TRUE)
writeRaster(ORN_rcl_buffer, 
				file.path(procDataDir, 
				paste0("ORN_reclass_", polygonBufferWidth, "km_buffered.tif")),
				overwrite = TRUE)
writeRaster(OHN_rcl_buffer, 
				file.path(procDataDir, 
				paste0("OHN_reclass_", polygonBufferWidth, "km_buffered.tif")),
				overwrite = TRUE)
writeRaster(urban_rcl_buffer, 
				file.path(procDataDir, 
				paste0("Urban_reclass_", polygonBufferWidth, "km_buffered.tif")),
				overwrite = TRUE)
writeRaster(rail_rcl_buffer, 
				file.path(procDataDir, 
				paste0("Rail_reclass_", polygonBufferWidth, "km_buffered.tif")),
				overwrite = TRUE)				
writeRaster(ORN_rclMajor, 
				file.path(procDataDir, "ORNMajor_reclass_FocalArea.tif"),
				overwrite = TRUE)	
writeRaster(ORN_rclMinor, 
				file.path(procDataDir, "ORNMinor_reclass_FocalArea.tif"),
				overwrite = TRUE)								

  # ESRI shapefile
st_write(ORN_buffer, 
				file.path(procDataDir, 
				paste0("ORN_", polygonBufferWidth, "km_buffered.shp")), 
				driver="ESRI Shapefile",
				append = FALSE)
st_write(OHN_buffer, 
				file.path(procDataDir, 
				paste0("OHN_", polygonBufferWidth, "km_buffered.shp")), 
				driver="ESRI Shapefile",
				append = FALSE)
st_write(urban_buffer, 
				file.path(procDataDir, 
				paste0("Urban_", polygonBufferWidth, "km_buffered.shp")), 
				driver="ESRI Shapefile",
				append = FALSE)
st_write(rail_buffer, 
				file.path(procDataDir, 
				paste0("Rail_", polygonBufferWidth, "km_buffered.shp")), 
				driver="ESRI Shapefile",
				append = FALSE)				

  # Save final output
  # Focal area
writeRaster(LULC_Focal, 
				file.path(procDataDir, "LULC_FocalArea.tif"),
				overwrite = TRUE)
writeRaster(LULC_vizFocal, 
				file.path(procDataDir, "LULC_Visualization_FocalArea.tif"),
				overwrite = TRUE)
  # Study area 							
writeRaster(LULC, 
				file.path(procDataDir, 
				paste0("LULC_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)
writeRaster(LULC_viz, 
				file.path(procDataDir,
				paste0("LULC_Visualization_", polygonBufferWidth, "km.tif")),
				overwrite = TRUE)

  # Buffered area
writeRaster(LULC_buffer, 
				file.path(procDataDir, 
				paste0("LULC_", polygonBufferWidth, "km_buffered.tif")), 
				overwrite=TRUE)
				
## End script--------------------------