#### a233: Royal Botanical Garden connectivity assessment
#### Script by Chloe Debyser

#### 3. Create resistance layer

############################################################################################################################
# This code:                                                                                                               #
# 1. Takes:                                                                                                                #
#    - The combined LULC layer produced in R code #2                                                                       #
#    - The LULC/resistance crosswalk                                                                                       #
# 2. Produces a resistance layer                                                                                           #
############################################################################################################################

#### Workspace ####
# Packages
library(tidyverse)
library(magrittr)
library(raster)

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

# Directories
setwd("E:/a233")
dataDir <- "Data/"

# Input parameters
bufferWidth <- 5 # In km

# Load data
      # Combined LULC layer
LULC <- raster(paste0(dataDir, "Processed/LULC_buffer", bufferWidth, "km.tif"))

      # Resistance crosswalk
crosswalk <- read.csv(paste0(dataDir, "Resistance crosswalks/ResistanceCrosswalk.csv"))

#### Resistance - Prepare resistance layer ####
# Reclassify
resistance <- LULC %>%
  reclassify(., rcl=crosswalk[, c("Destination_ID", "Resistance")])

# Save
writeRaster(resistance, paste0(dataDir, "Processed/Resistance_buffer", bufferWidth, "km.tif"))
