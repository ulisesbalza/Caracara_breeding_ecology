# Loadiing spatial data
rm(list=ls())

library(maptools)
library(spatstat)
library(rgdal)
library(sp)

data_fr <- read.csv("data/nidos_franklin_sitios.csv",header=T,sep=";",dec=".") #all sites 

# Window to spatial analysis 
S_fr <- readShapePoly("data/Penguin_poly_MCP_UTM.shp")
S_fr
SP_fr <- as(S_fr, "SpatialPolygons")
W_fr <- as(SP_fr, "owin")
plot(W_fr)
summary.owin(W_fr)

# Loading rockhopper penguin colonies
PPA_poly <- readShapePoly("data/PPA_poly.shp")
PPA <- as(PPA_poly, "SpatialPolygons")

#Check 
plot(W_fr, main = "Bahia Franklin")
plot(PPA, add=TRUE, col="red",border="red")

#Define data as PPP (planar point pattern) object

X_fr <- ppp(data_fr$x, data_fr$y,window=W_fr, unitname=c("metre","metres")) #all sites
#check
plot(X_fr, main = "Bahia Franklin, todos", pch = 19)
plot(PPA, add=TRUE, col="red", border="red")
summary(X_fr) #as it is all in meters, square units=square meters


