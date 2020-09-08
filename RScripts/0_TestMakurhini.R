#####################################################################
# a233 Test script - Makurhini library 
# 09-2020                                       					#
#                                                                   #                                       #
#                                                                   #
# Script by C Tucker for ApexRMS 									#
#####################################################################

# Workspace ---------------------------------------------------------

## Packages
  #To initially install
#require(devtools)
#install_github("connectscape/Makurhini", dependencies = TRUE, upgrade = "never")

library(Makurhini)

## Test of IIC and fractions measurements -----

  #Load data
data("vegetation_patches", package = "Makurhini")
nrow(vegetation_patches)
# 142

  #IIC
IIC <- MK_dPCIIC(nodes = vegetation_patches, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "IIC", distance_thresholds = 10000) #10 km
                
#Simple feature collection with 142 features and 5 fields
#geometry type:  POLYGON
#dimension:      XY
#bbox:           xmin: 3340120 ymin: 322869.6 xmax: 3739484 ymax: 696540.5
#CRS:            +proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs
# A tibble: 142 x 6
#      id   dIIC dIICintra dIICflux dIICconnector                                       geometry
# * <int>  <dbl>     <dbl>    <dbl>         <dbl>                                  <POLYGON [m]>
# 1     1 88.8    88.1      0.360           0.357 ((3676911 589967.3, 3676931 589895.5, 3676948…
# 2     2  0.736   0.0181   0.00766         0.710 ((3558044 696202.5, 3557972 696280.9, 3557957…
# 3     3  0.738   0.0119   0.0143          0.712 ((3569169 687776.4, 3569146 687749.5, 3569096…
# 4     4  0.719   0.00115  0.00194         0.716 ((3547316 685713.2, 3547362 685573.9, 3547390…
# 5     5  0.732   0.00554  0.0124          0.714 ((3567471 684357.4, 3567380 684214.3, 3567302…
# 6     6  0.732   0.0141   0.00677         0.711 ((3590569 672451.7, 3590090 672574.9, 3589912…
# 7     7  0.753   0.0352   0.0106          0.707 ((3570789 670959.4, 3570860 671015.4, 3570909…
# 8     8  0.724   0.00151  0.00642         0.716 ((3440118 666273.2, 3440372 666849.2, 3440584…
# 9     9  0.724   0.00163  0.00662         0.716 ((3451637 671232.4, 3451616 671287.1, 3451535…
# 10    10  0.732   0.00660  0.0111          0.714 ((3444396 671675.7, 3444715 671834.8, 3444873…
# … with 132 more rows

plot(IIC["dIIC"], breaks="jenks")


  # PC & fractions
PC <- MK_dPCIIC(nodes = vegetation_patches, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "PC", probability = 0.05,
                distance_thresholds = 10000)  
                
#PC
#               Simple feature collection with 142 features and 5 fields
#geometry type:  POLYGON
#dimension:      XY
#bbox:           xmin: 3340120 ymin: 322869.6 xmax: 3739484 ymax: 696540.5
#CRS:            +proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs
# A tibble: 142 x 6
#      id      dPC dPCintra dPCflux dPCconnector                                        geometry
# * <int>    <dbl>    <dbl>   <dbl>        <dbl>                                   <POLYGON [m]>
# 1     1 89.1     89.1     7.78e-4     0.       ((3676911 589967.3, 3676931 589895.5, 3676948 …
# 2     2  0.0194   0.0184  1.00e-3     5.72e-15 ((3558044 696202.5, 3557972 696280.9, 3557957 …
# 3     3  0.0152   0.0121  3.11e-3     3.82e-15 ((3569169 687776.4, 3569146 687749.5, 3569096 …
# 4     4  0.00153  0.00117 3.61e-4     5.05e-15 ((3547316 685713.2, 3547362 685573.9, 3547390 …
# 5     5  0.00833  0.00560 2.73e-3     0.       ((3567471 684357.4, 3567380 684214.3, 3567302 …
# 6     6  0.0143   0.0143  6.32e-5     0.       ((3590569 672451.7, 3590090 672574.9, 3589912 …
# 7     7  0.0358   0.0356  1.91e-4     1.57e-15 ((3570789 670959.4, 3570860 671015.4, 3570909 …
# 8     8  0.00228  0.00153 7.49e-4     2.89e-15 ((3440118 666273.2, 3440372 666849.2, 3440584 …
# 9     9  0.00216  0.00165 5.16e-4     5.17e-15 ((3451637 671232.4, 3451616 671287.1, 3451535 …
#10    10  0.00789  0.00668 1.21e-3     8.48e-15 ((3444396 671675.7, 3444715 671834.8, 3444873 …
# … with 132 more rows 

plot(PC["dPCflux"], breaks="jenks")

##--------------------------------------------------------------------------

## Test of IIC and fractions measurements for EMBL focal species-----
  # Additional packages
library(tidyverse)
library(raster)
library(sf)

 # Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A233 - Cootes to Escarpment"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")


## Run for EMBL, Blandings turtle------
species <- as.character("EMBL")

  # Load data
EMBL <- raster(file.path(outDir, paste0(species, "_Resistance_FocalArea.tif")))


  #IIC
IIC <- MK_dPCIIC(nodes = EMBL, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "IIC", distance_thresholds = 10000) #10 km
                plot(IIC)
#class      : RasterStack 
#dimensions : 750, 582, 436500, 5  (nrow, ncol, ncell, nlayers)
#resolution : 15, 15  (x, y)
#extent     : 1340475, 1349205, 11863365, 11874615  (xmin, xmax, ymin, ymax)
#crs        : +proj=lcc +lat_1=44.5 +lat_2=53.5 +lat_0=0 +lon_0=-85 +x_0=930000 +y_0=6430000 +ellps=GRS80 +units=m +no_defs 
#names      :           Id,         dIIC,    dIICintra,     dIICflux, dIICconnector 
#min values :  1.000000000,  0.998566458,  0.006963751,  0.991602708,   0.000000000 
#max values : 3.200000e+01, 8.358730e+01, 4.879439e+01, 3.479291e+01,  7.993606e-15 

plot(IIC[["dIIC"]]) #Test plot dIIC


  # PC & fractions
PC <- MK_dPCIIC(nodes = EMBL, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "PC", probability = 0.05,
                distance_thresholds = 10000)  
plot(PC[["dPCflux"]]) #Test plot dPCflux

#class      : RasterStack 
#dimensions : 750, 582, 436500, 5  (nrow, ncol, ncell, nlayers)
#resolution : 15, 15  (x, y)
#extent     : 1340475, 1349205, 11863365, 11874615  (xmin, xmax, ymin, ymax)
#crs        : +proj=lcc +lat_1=44.5 +lat_2=53.5 +lat_0=0 +lon_0=-85 +x_0=930000 +y_0=6430000 +ellps=GRS80 +units=m +no_defs 
#names      :          Id,         dPC,    dPCintra,     dPCflux, dPCconnector 
#min values : 1.000000000, 1.231267745, 0.006023568, 1.225244177,  0.000000000 
#max values :    32.00000,    85.82702,    42.20662,    43.62041,      0.00000 



## Run for BLBR------
species <- as.character("BLBR")

## Load data
BLBR <- raster(file.path(outDir, paste0(species, "_Resistance_FocalArea.tif")))


## Run for BLBR
  #IIC
IIC <- MK_dPCIIC(nodes = BLBR, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "IIC", distance_thresholds = 10000) #10 km
                plot(IIC)
plot(IIC[["dIIC"]]) #Test plot dIIC
                
#class      : RasterStack 
#dimensions : 750, 582, 436500, 5  (nrow, ncol, ncell, nlayers)
#resolution : 15, 15  (x, y)
#extent     : 1340475, 1349205, 11863365, 11874615  (xmin, xmax, ymin, ymax)
#crs        : +proj=lcc +lat_1=44.5 +lat_2=53.5 +lat_0=0 +lon_0=-85 +x_0=930000 +y_0=6430000 +ellps=GRS80 +units=m +no_defs 
#names      :           Id,         dIIC,    dIICintra,     dIICflux, dIICconnector 
#min values :  1.000000000,  0.942409948,  0.005587353,  0.936822595,   0.000000000 
#max values : 3.200000e+01, 5.411643e+01, 1.842403e+01, 3.569240e+01,  7.105427e-15 



  # PC & fractions
PC <- MK_dPCIIC(nodes = BLBR, attribute = NULL,
                distance = list(type = "centroid"),
                metric = "PC", probability = 0.05,
                distance_thresholds = 10000)  
plot(PC[["dPCflux"]]) #Test plot dPCflux

#class      : RasterStack 
#dimensions : 750, 582, 436500, 5  (nrow, ncol, ncell, nlayers)
#resolution : 15, 15  (x, y)
#extent     : 1340475, 1349205, 11863365, 11874615  (xmin, xmax, ymin, ymax)
#crs        : +proj=lcc +lat_1=44.5 +lat_2=53.5 +lat_0=0 +lon_0=-85 +x_0=930000 +y_0=6430000 +ellps=GRS80 +units=m +no_defs 
#names      :           Id,          dPC,     dPCintra,      dPCflux, dPCconnector 
#min values :  1.000000000,  0.754015632,  0.004329759,  0.749685874,  0.000000000 
#max values : 3.200000e+01, 5.834430e+01, 1.427717e+01, 4.406713e+01, 7.660539e-15 


