######################################################################
# a233 Royal Botanical Gardens Wildlife Corridor Mapping            
#  Headers, input parameter values   
#  09/2020                                                                          
#                                                                  
# Source file created by C Tucker  for ApexRMS 
#####################################################################


## Input parameters---------------------------------------------------------

  # Species
specieslist <- c("BLBR", "ODVI", "EMBL")

  # Final values used
  
polygonBufferWidth <- 20 # In km
# in km : the study area is made by buffering the polygon that encompases the Cootes to Escarpment EcoPark System

roadbuffer <- 12 #in m: Buffer road lines by constant so that rasterizing aligns with roads in SOLRIS

minCulvertSize <- 0.05
# About half of the OR values fall less than 0.015. According to Halton Road Ecology paper, values  must be > 0.05 for even smallest mammals    

culvertResistancevalue <- 10 # all suitable culverts given value of 10
bridgeResistancevalue <- 10 # all suitable bridges given value of 10

suitabilityThreshold <- 60
