######################################################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping                                             #
# Create study area from series of points and use it to create one LULC map                          #
#                                                                                                    #
# 1. Create study area                                                                               #
# - create polygon from points                                                                       #
# - Buffer polygon by 20km for study area                                                            #
# - Buffer study area by an additional 20% for connectivity analysis (to remove edge effects)        #
#                                                                                                    #
# 2. LULC crop and standardize rasters                                                               #
# - process all LULC layers: SOLRIS, OLCDB, ORN                                                      #
# - clip all LULC layers to buffered study area and reclassify                                       #
#                                                                                                    #
# 3. Combine all LULC rasters                                                                        #
#                                                                                                    #
# 4. Clip all intermediate and final rasters to study area for visualizations                        #
#                                                                                                    #
# Script created by Bronwyn Rayfield and Chloé Debeyser for ApexRMS                                  #
######################################################################################################

# Workspace ------------------------
# Packages
library(tidyverse)
library(sf)
library(raster)

# Directories
projectDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A233_RBGConnectivity/a233"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")

# Input parameters
polygonBufferWidth <- 20 # in km : the study area is made by buffering the polygon that encompases the Cootes to Escarpment EcoPark System

# Read in data
  # Study area
studyAreaPoints <- read_csv(file.path(paste0(dataDir, "/Study Area"), "Vertices 20200610 DAG.csv")) # Study area vertices

  # LULC data
SOLRISRaw <- raster(file.path(paste0(dataDir, "/Land use land cover/SOLRIS_Version_3_0"), "SOLRIS_Version_3_0_LAMBERT.tif"))
OLCDBRaw <- raster(file.path(paste0(dataDir, "/Land use land cover/Provincial Land Cover 1996 - 28 Class Grid"), "plc1996_28.tif"))
ORNRaw <- st_read(file.path(paste0(dataDir, "/Land use land cover/Ontario_Road_Network__ORN__Segment_With_Address-shp"), "Ontario_Road_Network__ORN__Segment_With_Address.shp"))
OHNRaw <- st_read(file.path(paste0(dataDir, "/Land use land cover/Ontario_Hydro_Network__OHN__-_Waterbody-shp"), "Ontario_Hydro_Network__OHN__-_Waterbody.shp"))
urbanRaw <- st_read(file.path(paste0(dataDir, "/Land use land cover/Built-Up_Area-shp"), "Built-Up_Area.shp"))
  # Resistance crosswalk
crosswalk <- read_csv(file.path(paste0(dataDir, "/Resistance"), "GenericSpeciesResistanceCrosswalk.csv"))

# Create study area -------------------------------------
# Create spatial polygon from points
polygon <- studyAreaPoints %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# Reproject polygon to the same crs as lulcRaw
polygonProjected <- st_transform(polygon, crs(SOLRISRaw))
# polygonProjected <- st_read(paste(studyAreaDir, "polygon_projected.shp", sep="/"))

# Create buffer around polygon for the study area using polygon buffer width
studyArea <- st_buffer(polygonProjected, (polygonBufferWidth*1000))

# Create bounding box around study area
# Caculate 20% of the average of the max length and width of the bounding box to be the study area buffer width
# Connectivity analyses will be run in this buffered study area to minimize edge effects
studyAreaBBox <- st_bbox(studyArea)
studyAreaBufferWidth <- round(((studyAreaBBox$ymax - studyAreaBBox$ymin) + (studyAreaBBox$xmax - studyAreaBBox$xmin))/2 * 0.2, digits=0) #in m
studyAreaBuffer <- st_buffer(studyArea, studyAreaBufferWidth)


# LULC - Crop data to STUDY AREA BUFFER and standardize rasters ---------------------
# SOLRIS 
# Crop to buffered study area
SOLRIS_buffer <- SOLRISRaw %>%
  crop(., extent(studyAreaBuffer), snap="out") %>% # Crop SOLRIS to study area buffer extent
  mask(., mask=studyAreaBuffer) %>% # Clip to buffered study area
  trim(.) # Trim extra white spaces
# unique(SOLRIS)

# Reclassify
SOLRIS_rcl_buffer <- SOLRIS_buffer %>%
  reclassify(., rcl=crosswalk[which(crosswalk$Source_Name == "SOLRIS"), c("Source_ID", "Destination_ID")])

# OLCDB
# Align with SOLRIS
# Buffer SOLRIS extent
ext_OLCDB <- extend(extent(SOLRIS_buffer), 1000)

# Crop to buffered study area
OLCDB_buffer <- OLCDBRaw %>% 
  crop(., y=ext_OLCDB, snap="out") %>% # Crop OLCDB to new extent 
  projectRaster(from=., crs=crs(SOLRIS_buffer), method='ngb') %>% # Project to SOLRIS CRS 
  resample(., SOLRIS_buffer, method='ngb') %>% # Resample to SOLRIS resolution
  mask(., mask=studyAreaBuffer) %>% # Clip to buffered study area
  trim(.) # Trim extra white spaces
# unique(OLCDB)

# Reclassify
OLCDB_rcl_buffer <- OLCDB_buffer %>%
  reclassify(., rcl=crosswalk[which(crosswalk$Source_Name == "OLCDB"), c("Source_ID", "Destination_ID")])

# ORN
# Crop to buffered study area
ORN_buffer <- ORNRaw %>%
  st_transform(., crs=st_crs(SOLRIS_buffer)) %>% # Transform to SOLRIS CRS
  st_intersection(., studyAreaBuffer) %>% # Clip to buffered study area
  st_buffer(., 12) # Buffer road lines by 12m so that rasterizing aligns with roads in SOLRIS

# Reclassify
ORN_rcl_buffer <- ORN_buffer %>%
  left_join(., crosswalk[which(crosswalk$Source_Name == "ORN"), c("Destination_ID", "Source_Label")], by=c("ROAD_CLASS" = "Source_Label")) %>% # Add Destination_ID
  rasterize(., SOLRIS_buffer, field="Destination_ID") # Rasterize

# OHN
# Crop to buffered study area
OHN_buffer <- OHNRaw %>%
  st_transform(., crs=st_crs(SOLRIS_buffer)) %>% # Transform to SOLRIS CRS
  mutate(AREA = as.numeric(st_area(.))) %>% # Calculate the area of each waterbody
  dplyr::filter(AREA >= 100) %>% # Keep waterbodies larger than 100m2
  st_intersection(., studyAreaBuffer) # Clip to buffered study area

OHN_rcl_buffer <- OHN_buffer %>%
  left_join(., crosswalk[which(crosswalk$Source_Name == "OHN"), c("Destination_ID", "Source_Label")], by=c("WATERBODY_" = "Source_Label")) %>% # Add Destination_ID
  rasterize(., SOLRIS_buffer, field="Destination_ID") # Rasterize

# Urban Areas
urban_buffer <- urbanRaw %>%
  st_transform(., crs=st_crs(SOLRIS_buffer)) %>% # Transform to SOLRIS CRS
  st_intersection(., studyAreaBuffer)

urban_rcl_buffer <- urban_buffer %>%
  #  filter(COMMUNIT_1 == "Built-up Area Impervious") %>%
  left_join(., crosswalk[which(crosswalk$Source_Name == "Built-up Areas"), c("Destination_ID", "Source_Label")], by=c("COMMUNIT_1" = "Source_Label")) %>% # Add Destination_ID
  rasterize(., SOLRIS_buffer, field="Destination_ID") # Rasterize


# LULC - Combine rasters ---------------------------
# Remove undifferentiated cells from SOLRIS so that they can take the value from OLCDB
SOLRIS_rcl_buffer[SOLRIS_rcl_buffer == 250] <- NA

# Merge rasters: urban gets priority, followed by ORN, OHN, SOLRIS, and last OLCDB
LULC_buffer <- merge(urban_rcl_buffer, ORN_rcl_buffer, OHN_rcl_buffer, SOLRIS_rcl_buffer, OLCDB_rcl_buffer)


# Crop data to STUDY AREA ------------------------------
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

# Urban
urban <- urban_buffer %>%
  st_intersection(., studyArea) # Clip to study area

urban_rcl <- urban_rcl_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop SOLRIS to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.) # Trim extra white spaces

# LULC
LULC <- LULC_buffer %>%
  crop(., extent(studyArea), snap="out") %>% # Crop to study area extent
  mask(., mask=studyArea) %>% # Clip to study area
  trim(.)

# Save outputs --------------------------
# Save intermediate output
st_write(polygonProjected, file.path(outDir, "FocalArea.shp"), driver="ESRI Shapefile")
st_write(studyArea, file.path(outDir, paste0("StudyArea_", polygonBufferWidth, "km.shp")), driver="ESRI Shapefile")
st_write(studyAreaBuffer, file.path(outDir, paste0("StudyArea_", polygonBufferWidth, "km_buffered.shp")), driver="ESRI Shapefile")
# Study Area
writeRaster(SOLRIS, file.path(outDir, paste0("SOLRIS_", polygonBufferWidth, "km.tif")))
writeRaster(SOLRIS_rcl, file.path(outDir, paste0("SOLRIS_reclass_", polygonBufferWidth, "km.tif")))
writeRaster(OLCDB, file.path(outDir, paste0("OLCDB_", polygonBufferWidth, "km.tif")))
writeRaster(OLCDB_rcl, file.path(outDir, paste0("OLCDB_reclass_", polygonBufferWidth, "km.tif")))
st_write(ORN, file.path(outDir, paste0("ORN_", polygonBufferWidth, "km.shp")), driver="ESRI Shapefile")
writeRaster(ORN_rcl, file.path(outDir, paste0("ORN_reclass_", polygonBufferWidth, "km.tif")))
st_write(OHN, file.path(outDir, paste0("OHN_", polygonBufferWidth, "km.shp")), driver="ESRI Shapefile")
writeRaster(OHN_rcl, file.path(outDir, paste0("OHN_reclass_", polygonBufferWidth, "km.tif")))
st_write(urban, file.path(outDir, paste0("Urban_", polygonBufferWidth, "km.shp")), driver="ESRI Shapefile")
writeRaster(urban_rcl, file.path(outDir, paste0("Urban_reclass_", polygonBufferWidth, "km.tif")))
# Study Area Buffer
writeRaster(SOLRIS_buffer, file.path(outDir, paste0("SOLRIS_", polygonBufferWidth, "km_buffered.tif")))
writeRaster(SOLRIS_rcl_buffer, file.path(outDir, paste0("SOLRIS_reclass_", polygonBufferWidth, "km_buffered.tif")))
writeRaster(OLCDB_buffer, file.path(outDir, paste0("OLCDB_", polygonBufferWidth, "km_buffered.tif")))
writeRaster(OLCDB_rcl_buffer, file.path(outDir, paste0("OLCDB_reclass_", polygonBufferWidth, "km_buffered.tif")))
st_write(ORN_buffer, file.path(outDir, paste0("ORN_", polygonBufferWidth, "km_buffered.shp")), driver="ESRI Shapefile")
writeRaster(ORN_rcl_buffer, file.path(outDir, paste0("ORN_reclass_", polygonBufferWidth, "km_buffered.tif")))
st_write(OHN_buffer, file.path(outDir, paste0("OHN_", polygonBufferWidth, "km_buffered.shp")), driver="ESRI Shapefile")
writeRaster(OHN_rcl_buffer, file.path(outDir, paste0("OHN_reclass_", polygonBufferWidth, "km_buffered.tif")))
st_write(urban_buffer, file.path(outDir, paste0("Urban_", polygonBufferWidth, "km_buffered.shp")), driver="ESRI Shapefile")
writeRaster(urban_rcl_buffer, file.path(outDir, paste0("Urban_reclass_", polygonBufferWidth, "km_buffered.tif")))

# Save final output
writeRaster(LULC, file.path(outDir, paste0("LULC_", polygonBufferWidth, "km.tif")))
writeRaster(LULC_buffer, file.path(outDir, paste0("LULC_", polygonBufferWidth, "km_buffered.tif")), overwrite=TRUE)

