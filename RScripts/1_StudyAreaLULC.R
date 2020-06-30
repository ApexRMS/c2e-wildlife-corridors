######################################################################################################
# Create study area from series of points and use it to create one LULC map                          #
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
# Script created by Bronwyn Rayfield and Chloé Debeyser for ApexRMS                                  #
######################################################################################################

# Workspace ------------------------
# Packages
library(tidyverse)
library(sf)
library(raster)

# Directories
projectDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A233_RBGConnectivity"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")

# Input parameters
polygonBufferWidth <- 20 # in km : the study area is made by buffering the polygon that encompases the Cootes to Escarpment EcoPark System

# Read in data
  # Study area
studyAreaPoints <- read_csv(file.path(paste0(dataDir, "/Study Area"), "Vertices 20200610 DAG.csv")) # Study area vertices

  # LULC data
OLCDBRaw <- raster(file.path(paste0(dataDir, "/Land use land cover/Provincial Land Cover 1996 - 28 Class Grid"), "plc1996_28.tif"))
ORNRaw <- st_read(file.path(paste0(dataDir, "/Land use land cover/Ontario_Road_Network__ORN__Segment_With_Address-shp"), "Ontario_Road_Network__ORN__Segment_With_Address.shp"))
SOLRISRaw <- raster(file.path(paste0(dataDir, "/Land use land cover/SOLRIS_Version_3_0"), "SOLRIS_Version_3_0_LAMBERT.tif"))

  # Resistance crosswalk
crosswalk <- read.csv(file.path(paste0(dataDir, "/Resistance"), "ResistanceCrosswalk.csv"))

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


# LULC - Crop and standardize rasters ---------------------
# SOLRIS
# Crop
SOLRIS <- SOLRISRaw %>%
  crop(., extent(studyAreaBuffer), snap="out") %>% # Crop SOLRIS to study area buffer extent
  mask(., mask=studyAreaBuffer) %>% # Clip to buffered study area
  trim(.) # Trim extra white spaces
# unique(SOLRIS)

# Reclassify
SOLRIS_rcl <- SOLRIS %>%
  reclassify(., rcl=crosswalk[which(crosswalk$Source_Name == "SOLRIS"), c("Source_ID", "Destination_ID")])

# OLCDB
# Align with SOLRIS
# Buffer SOLRIS extent
ext_OLCDB <- extend(extent(SOLRIS), 1000)

# Crop
OLCDB <- OLCDBRaw %>% 
  crop(., y=ext_OLCDB, snap="out") %>% # Crop OLCDB to new extent 
  projectRaster(from=., crs=crs(SOLRIS), method='ngb') %>% # Project to SOLRIS CRS 
  resample(., SOLRIS, method='ngb') %>% # Resample to SOLRIS resolution
  mask(., mask=studyAreaBuffer) %>% # Clip to buffered study area
  trim(.) # Trim extra white spaces
# unique(OLCDB)

# Reclassify
OLCDB_rcl <- OLCDB %>%
  reclassify(., rcl=crosswalk[which(crosswalk$Source_Name == "OLCDB"), c("Source_ID", "Destination_ID")])

# ORN
# Crop
ORN <- ORNRaw %>%
  st_transform(., crs=st_crs(SOLRIS)) %>% # Transform to SOLRIS CRS
  st_intersection(., studyAreaBuffer) # Clip to buffered study area

# Reclassify
ORN_rcl <- ORN %>%
  left_join(., crosswalk[which(crosswalk$Source_Name == "ORN"), c("Destination_ID", "Source_Label")], by=c("ROAD_CLASS" = "Source_Label")) %>% # Add Destination_ID
  rasterize(., SOLRIS, field="Destination_ID") # Rasterize

# LULC - Combine rasters ---------------------------
# Remove undifferentiated cells from SOLRIS so that they can take the value from OLCDB
SOLRIS_rcl[SOLRIS_rcl == 250] <- NA

# Merge rasters: ORN gets priority, followed by SOLRIS, followed by OLCDB
LULC <- merge(ORN_rcl, SOLRIS_rcl, OLCDB_rcl)

# Save outputs --------------------------
# Save intermediate output
st_write(polygonProjected, file.path(outDir, "polygon_projected.shp"), driver="ESRI Shapefile")
st_write(studyArea, file.path(outDir, "studyarea_unbuffered.shp"), driver="ESRI Shapefile")
st_write(studyAreaBuffer, file.path(outDir, "studyarea_buffered.shp"), driver="ESRI Shapefile")
writeRaster(SOLRIS, file.path(outDir, paste0("SOLRIS_", polygonBufferWidth, "km_buffered.tif")))
writeRaster(SOLRIS_rcl, file.path(outDir, paste0("SOLRIS_reclass_", polygonBufferWidth, "km_buffered.tif")))
writeRaster(OLCDB, file.path(outDir, paste0("OLCDB_", polygonBufferWidth, "km_buffered.tif")))
writeRaster(OLCDB_rcl, file.path(outDir, paste0("OLCDB_reclass_", polygonBufferWidth, "km_buffered.tif")))
st_write(ORN, file.path(outDir, paste0("ORN_", polygonBufferWidth, "km_buffered.shp")), driver="ESRI Shapefile")
writeRaster(ORN_rcl, file.path(outDir, paste0("ORN_reclass_", polygonBufferWidth, "km_buffered.tif")))

# Save final output
writeRaster(LULC, file.path(outDir, paste0("LULC_", polygonBufferWidth, "km_buffered.tif")))

