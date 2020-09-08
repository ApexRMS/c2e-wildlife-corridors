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
