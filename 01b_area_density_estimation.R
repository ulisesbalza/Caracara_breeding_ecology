rm(list=ls())

#AREA TO WHICH REFER THE DENSITY ESTIMATION, BASED ON DISTANCE TO THE ROCKHOPPER PENGUIN COLONY

library(maptools)
library(spatstat)
library(rgdal)
library(raster)
library(sp)
library(geosphere)
library(rgeos)
library(sf)
library(RColorBrewer)
library(ggplot2)

# Loading data
data_fr_2018 <- read.csv("data/nidos_franklin_2018.csv",header=T,sep=";",dec=".")

# Spatial window
S_fr <- readShapePoly("data/Penguin_poly_MCP_UTM.shp")
S_fr
SP_fr <- as(S_fr, "SpatialPolygons")
W_fr <- as(SP_fr, "owin")
plot(W_fr)
summary.owin(W_fr)

# Rockhopper penguin colonies
PPA_poly <- readShapePoly("data/PPA_poly.shp")
PPA <- as(PPA_poly, "SpatialPolygons")

# Point pattern of the data
X_fr_2018 <- ppp(data_fr_2018$x, data_fr_2018$y,window=W_fr, unitname=c("metre","metres"))

# line transects
x <- c(393729.32, 393138.26, 392248.93, 392623.99)
y <- c(3918146.05-100,3918101.66-100, 3919743-100, 3920058.02-100)
# make SpatialPoints
points <- sp::SpatialPoints(cbind(x,y))
# use as to convert to line
sp_line <- as(points,"SpatialLines")


# Check everything
plot(X_fr_2018, main = "Bahia Franklin 2018", pch = 19, cex=2)
plot(PPA, add=T, col="red", border="red")
plot(sp_line, add=T, lwd=10, col="lightblue")


# Calculate distance from transect to the colony

#Generate regular points on the transect
ptsreg <- spsample(sp_line, 30, type = "regular")

# Create a raster of distance to the rockhopper colony
#create a fishnet Cellsize in m
# First need both files in sf
franklinsf <- (st_as_sf(S_fr))
ppasf <- (st_as_sf(PPA))

# Make the grid
grid <- st_make_grid(franklinsf, cellsize = 105, what = "centers") #cellsize= effective band size (from Distance)

plot(grid, pch=20, main="Grid every 105 m")
plot(ppasf, add=T, col="red", borders="red")

# calculate distances from points to PPA
#transform from polygon shape to line
ppasf <- st_cast(ppasf, "MULTILINESTRING")

#Calculation of distance to nearest rockhopper patch
ppa_sp <- as(ppasf, 'Spatial')
grid_sp <- as(grid, 'Spatial')

dist <- apply(gDistance(grid_sp, ppa_sp, byid=TRUE),2,min)

#plotting
#create a data.frame with the distance and the coordinates of the points
df <- data.frame(dist = as.vector(dist),
                 st_coordinates(grid))

col_dist <- brewer.pal(11, "RdGy")

ggplot(df, aes(X, Y, fill = dist))+ #variables
  geom_tile()+ #geometry
  scale_fill_gradientn(colours = rev(col_dist))+ #colors for plotting the distance
  labs(fill = "Distance (m)")+ #legend name
  theme_void()+ #map theme
  theme(legend.position = "bottom") #legend position


coordinates(df) <- ~ X + Y
# coerce to SpatialPixelsDataFrame
gridded(df) <- TRUE
# coerce to raster
ppa_raster <- raster(df)

plot(ppa_raster, main="Distance to rockhopper (m)")
plot(PPA, add=T)
plot(W_fr, add=T, lwd=2)
plot(sp_line, add=T, lwd=10, col="lightblue")
plot(ptsreg, pch=10, add=T)

# Extraer raster data
buffer_105 <- extract(ppa_raster, ptsreg, buffer=105, fun=mean, df=TRUE)
summary(buffer_105$dist) 

# what amount of area corresponds with this maximum distances?

m <- c(0, 288, 1)
m <- matrix(m, ncol=3, byrow = T)

ppa_raster_filtrado <- reclassify(ppa_raster, m, right = T)
plot(ppa_raster, main="Distance to rockhopper (m)")
plot(ppa_raster_filtrado)
plot(W_fr, add=T, lwd=2)
plot(ptsreg, pch=10, add=T)

#clip raster
franklin <- as(W_fr, "SpatialPolygons")
clip_raster <- crop(ppa_raster_filtrado, extent(franklin))
ppa_raster_filtrado_cortado <- mask(clip_raster, franklin)

plot(ppa_raster_filtrado_cortado)
plot(franklin, add=T, lwd=2)

# count cells =1 inside the window

head(freq(ppa_raster_filtrado_cortado)) #366 celdas
(366 *(105*105))/10000 #total area in hectares

area(franklin)/10000 #window area in hectares

403.515/587.5996 # relative area to which refer the density estimation
