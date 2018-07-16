#####################################
## Interpolate E-OBS mean daily average sea level pressure 
## (units = hecto Pascals) to Irish hectad scale
## 
## inputs:  * eobs.RData - data prepared by Jon Yearsley and shared in the 
##            Grassland dropbox folder, summer 2017
##          * ireland_coastline shapefile (I got this from Jon Yearsley but I
##                think it might be a naturalEarth file originally)
##          * IE_10km_hecs.shp shapefile of Irish hectads created in QGIS by wg
## outputs: krg_mean_pp - a spatialGridDataFrame holding 
##          interpolated values of 98th quantile of maximum summer temperature 
##          over the years 1995-2016 at the hectad scale.  
##          This includes many ocean blocks where 
##          predictions are likely poor.  This should be masked before being
##          used for analysis.  This is the object to use for most things.
##          
##          krg_pp_rast - a rasterized version of krg_mean_pp  
##          
##          krig_pp_map - a rasterBrick version of krg_mean_pp,
##          useful for easy plotting (plot method shows variance by default). 
##          krg_summer_tx_predict is probably the object to use for analysis.
##
## TODO:  - Perhaps try the pseudo-matlab version Jon mentioned?
##          alternatviely, try approxfun() or the R equivalent of Matlab 
##          "interp" "interp2" for 2d.
## 
## author: Willson Gaul
## created: 13 Dec 2017
## last modified: 3 May 2018
#############################################


setwd("~/Documents/Data_Analysis/UCD/predictor_variables/eobs")
library(wgutil)
library(Hmisc)
library(rgdal)
library(gdalUtils)
library(ggplot2)
library(fields)
library(gstat)
library(raster)
library(tidyverse)
load("../../data/eobs.RData")

## ----------------------- load coastline for masking -------------------------
# Load Ireland coastline
ir <- readOGR(dsn='../../mapping/data/', layer='ireland_coastline')
ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))

### ----------------- prepare hectad raster -----------------------------------
# rasterize the hectad shapefile so that I can predict to it later
src_filename <- "../../qgis/find_IE_hectads/data/IE_10km_hecs.shp"
dst_filename <- "IE_10km_hecs_raster.tif"
irish_hec_raster <- gdalUtils::gdal_rasterize(src_filename, 
                                              dst_filename = dst_filename, 
                                              a = "hctd_cd", 
                                              tr = c(10000, 10000),
                                              output_Raster = F)
irish_hec_raster <- raster::raster("IE_10km_hecs_raster.tif")
irish_hec_raster <- projectRaster(from = irish_hec_raster, 
                                  crs = CRS("+init=epsg:29903"))
### --------------------- end hectad raster -----------------------------------

## calculate mean pp over all years -------------------------------------------
# tx is maximum daily temperature (deg C). 
eobs_data$location <- paste(eobs_data$latitude, eobs_data$longitude, sep = "_")
mean_pp <-  select(eobs_data, longitude, latitude, 
                             location, pp) %>%
  group_by(location) %>%
  mutate(mean_pp = mean(pp)) %>%
  select(-pp) %>%
  distinct()
## end calculating pp -------------------------------------

## make into spatial df ------------------------------------------------------
spat_mean_pp <- mean_pp[, 3:ncol(mean_pp)]
coordinates(spat_mean_pp) <- as.matrix(mean_pp[, 1:2])
proj4string(spat_mean_pp) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_mean_pp <- spTransform(spat_mean_pp, 
                            CRS("+init=epsg:29903"))
dimnames(spat_mean_pp@coords)[[2]][which(
  dimnames(spat_mean_pp@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_mean_pp@coords)[[2]][which(
  dimnames(spat_mean_pp@coords)[[2]] == "latitude")] <- "northings"

## ----------------------- interpolate by krigging ----------------------------
# Interpolate mean_pp to the raster grid irish_hec_raster
irish_spat_grid <- as(irish_hec_raster, 'SpatialGrid')

# create empirical variogram for mean pp
gs_mean_pp <- gstat(id = "mean_pp", 
                    formula = mean_pp ~ 1, 
                    data = spat_mean_pp) # make gstat object
v_mean_pp <- variogram(gs_mean_pp)
#vgm(as.character(vgm()[, 1])) # show vgm options
f_var_mean_pp <- fit.variogram(v_mean_pp, 
                               vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", 
                                     "Ste"))) 
f_var_mean_pp
plot(variogramLine(f_var_mean_pp, max(v_mean_pp$dist)), type = 'l')
points(v_mean_pp[, 2:3], pch = 20, col = 'red')

# use variogram in kriging interpolation
krg_mean_pp <- gstat(formula = mean_pp~1, 
                     data = spat_mean_pp, model = f_var_mean_pp)

krg_mean_pp <- predict(krg_mean_pp, irish_spat_grid)
names(krg_mean_pp) <- c("mean_pp", "variance")

# make rasters
krg_mean_pp_rast <- raster::raster(krg_mean_pp)

# masking with shapefile, but I think the shapefile cuts off some hectads that
# do have data (some sp. datasets have 1014 hectads)
plot(raster::mask(krg_mean_pp_rast, ir_TM75))

krig_mean_pp_map <- raster::brick(krg_mean_pp)
krig_mean_pp_map <- raster::mask(krig_mean_pp_map, ir_TM75) 
names(krig_mean_pp_map) <- c("mean_pp", "variance")
plot(krig_mean_pp_map, main = "Mean pp")


### save outputs ---------------------------------------------------------------
save(krg_mean_pp, krg_mean_pp_rast, krig_mean_pp_map, 
     file = "mean_pp_hectad.RData")