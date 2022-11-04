rm(list=ls())

library(rgeos)
library(tidyverse)
library(sf)
library(raster)
library(spdep)
library(cubature)
library(MuMIn)

#---------------------------------#
# Descriptive plots ####
#---------------------------------#

plot(density(X_fr), main= "Kernel density", riblab="Density (nests/km2)", 
     ribscale=1*1000000, ribside="left", ribn=5)
plot(PPA_poly, add=T, col="black")
plot(X_fr, add=T, pch = 19, cex=2)
plot(W_fr, add=T, lwd=2)

# Nearest neighbors
d_fr <- nndist(X_fr)
summary(d_fr) 

# First-order density analysis

# First, we are going to create a point process pattern based on distance to the penguin colony ####--
# Create a raster of distances to the colony

#create a fishnet Cellsize in m
franklinsf <- (st_as_sf(S_fr))
ppasf <- (st_as_sf(PPA))

#make grid
grid <- st_make_grid(franklinsf, cellsize = 30, what = "centers")

#check
plot(grid, pch=20, main="Grid 30 m")
plot(ppasf, add=T, col="red", borders="red")


# calculate distances from points to penguin colony
# transform from polygon shape to line
ppasf <- st_cast(ppasf, "MULTILINESTRING")

#Calculation of distance to nearest rockhopper patch
ppa_sp <- as(ppasf, 'Spatial')
grid_sp <- as(grid, 'Spatial')

dist <- apply(gDistance(grid_sp, ppa_sp, byid=TRUE),2,min)

#create a data.frame with the distance and the coordinates of the points
df <- data.frame(dist = as.vector(dist),
                 st_coordinates(grid))
coordinates(df) <- ~ X + Y
# coerce to SpatialPixelsDataFrame
gridded(df) <- TRUE
# coerce to raster
ppa_raster <- raster(df)
ppa_raster <- mask(ppa_raster, SP_fr)
plot(ppa_raster, main="Distance to rockhopper (m)")
plot(PPA, add=T)
plot(W_fr, add=T, lwd=2)

#raster into an image covariate that spatstat can read
distppa_im <- as.im(ppa_raster)


# Loading topographic variables (elevation, aspect and slope) ####

rr_slope<- raster("data/franklin_pendiente.tif")
raster_slope <- as.im(rr_slope)

rr_aspect<- raster("data/franklin_orientacion.tif")
raster_aspect <- as.im(rr_aspect)

rr_elevation<- raster("data/franklin_elevacion.tif")
raster_elevation <- as.im(rr_elevation)


# Set the candidate models
# Model list

# Linear models
pp.int <- ppm(X_fr, ~ 1) #no trend (homogeneous)
pp.x <- ppm(X_fr, ~ x) #linear with lon
pp.y <- ppm(X_fr, ~ y) #linear with lat
pp.xy <- ppm(X_fr, ~ x + y) #linear trend lat + lon
pp.rockhoper <- ppm(X_fr, ~ distppa_im) #linear with distance to the colony
# GAM models (optimal may be intermediate)
pp.slope <- ppm(X_fr, ~ s(raster_slope), use.gam = T)
pp.aspect <- ppm(X_fr, ~ s(raster_aspect), use.gam = T)
pp.elevation <- ppm(X_fr, ~ s(raster_elevation), use.gam = T)


# Model selection by AICc
# AICc of models
models <-  data.frame(model = c("null", "lon", "lat", "lonlat", 
                     "rockhopper", "slope", "aspect", "elevation"),
             AICc = c(AICc(pp.int), AICc(pp.x), AICc(pp.y), AICc(pp.xy), 
                      AICc(pp.rockhoper), AICc(pp.slope), AICc(pp.aspect), AICc(pp.elevation)))

models$deltaAICc <- models$AICc-620.93 #AICc of null model
head(models)

summary_models <- models[order(models$AICc),]
summary_models



#Best model summary and diagnostics
summary(pp.rockhoper) 
diagnose.ppm(pp.rockhoper)


# Spatial prediction based on the best model
rockhopper_predict <-  predict(pp.rockhoper, type = "trend")

plot(rockhopper_predict, riblab="Density (nests/km2)", 
     ribscale=1*1000000, ribside="left", col=rev(heat.colors(50)), ribn=5, main="Distance to rockhopper colony model")
plot(W_fr, add=T, lwd=2)
contour(rockhopper_predict, add = TRUE, levels = c(0.0000039), labels = "", #show densities higher than the mean
        labcex = 1) #marco la media
plot(X_fr, add=T, pch=19, cex=1.5)
plot(PPA, add=T)


# spatial simulation with best model
simu_rockhopper <- rmpoispp(rockhopper_predict, types = c("caracara"), nsim=4, drop=F)
plot(simu_rockhopper, pch=19, cex=2, main = "Simulation test")



# K statistics based on heterogeneous patterns
# We'll use the best first-order density to analyze if there is evidence of clustering
# or repulsion AFTER taking into account variation in density driven by rockhopper penguins
# (i.e., something that it is not related to interaction among events-nests)

# pattern vs simulated as expected
L_ppa_csr <- envelope(X_fr, Linhom, nsim = 99, rank = 2,
                   correction = "none", simulate =
                     expression(rpoispp(rockhopper_predict)), global = F)

plot(L_ppa_csr, main="Franklin Bay, inhomogeneous pattern by rockhopper distance", legend = F)

# g statistics inhom

P_ppa_csr <- envelope(X_fr, pcf, nsim = 99, rank = 2,
                      correction = "none", simulate =
                        expression(rpoispp(rockhopper_predict)), global = F)
plot(P_ppa_csr, main="Franklin Bay, inhomogeneous pattern by rockhopper distance", legend = F)


#NN distances inhomo
G_ppa_env_fr <- envelope(X_fr, Gest, nsim = 99, rank = 2,
                         correction = "none", simulate =
                           expression(rpoispp(rockhopper_predict)), global = F)
plot(G_ppa_env_fr, main="G(nn) inhomogeneous pattern by rockhopper distance", legend = FALSE)



# Homogeneous patterns

#L(r)
plot(env_L_fr, main="L(r) Homogeneous", legend = F) 

# g
plot(Penv_fr, ylim=c(0,20), main="g(r) Homogeneous", legend = FALSE)

#G_NN
Genv_fr <- envelope(X_fr, Gest, nsim = 99, rank = 2,
                    correction = "none", global = FALSE)
plot(Genv_fr, main="G(nn) Homogeneous", legend = FALSE)

