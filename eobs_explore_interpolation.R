###################################
## Interpolate E-OBS climate data (1995-2016) to Irish hectad scale
## 
## This is an exploratory script, trying to figure out how to summarize 
## the E-OBS data at different temporal scales before interpolating.
## 
## inputs: eobs.Rdata - data prepared by Jon Yearsley and shared in the 
##            Grassland dropbox folder, summer 2017
## outputs: individual climate variables at the location and scale of 
##            Irish hectads
## 
## author: Willson Gaul
## created: 25 Sep 2017
## last modified: 2 Nov 2017
#################################### 

setwd("~/Documents/Data_Analysis/UCD/predictor_variables/eobs")
library(wgutil)
library(rgdal)
library(gdalUtils)
library(ggplot2)
library(fields)
library(gstat)
library(raster)
rm(list = ls())
load("../../data/eobs.RData")

## Status
# This now can produce the result format I want.  Next step is to figure out exactly what I want.  Interpolate for ever day so I can match observations to variable values?  Average values by week or month?

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

## ------------------ convert eobs lat/long to OSI east/north ----------------
# make eobs_data a spatial data frame
spat_eobs <- eobs_data[, 3:14]
coordinates(spat_eobs) <- as.matrix(eobs_data[, 1:2])
# set proj4string based on WGS84 here: http://spatialreference.org/ref/epsg/4326/
proj4string(spat_eobs) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_eobs <- spTransform(spat_eobs, CRS("+init=epsg:29903"))

## -------------------------- end coordinate conversion ------------------------

## ------------------------- interpolate --------------------------------------
# I want to interpolate data in spat_eobs to the raster grid irish_hec_raster
# following this: http://rspatial.org/analysis/rst/4-interpolation.html
irish_spat_grid <- as(irish_hec_raster, 'SpatialGrid')

# create empirical variogram for jan 1 2015
jan1 <- spat_eobs[spat_eobs$date == "2015-01-01", ]
gs_eobs <- gstat(id = "jan1_2015_temp", formula = tg~1, 
                 data = jan1) # make gstat object
v <- variogram(gs_eobs, width = 1000)

f_var <- fit.variogram(v, vgm("Exp"))
f_var
plot(variogramLine(f_var, max(v$dist)), type = 'l', ylim = c(0, 0.5))
points(v[, 2:3], pch = 20, col = 'red')

# use variogram in kriging interpolation
krg <- gstat(formula = tg~1, data = jan1, model = f_var)

krg_predict <- predict(krg, irish_spat_grid)
#spplot(krg_predict)

## krg_predict has the results I want
krg_rast <- raster::raster(krg_predict)

# masking with shapefile, but I think the shapefile cuts off some hectads that
# do have data (some sp. datasets have 1014 hectads)
results_masked <- raster::mask(krg_rast, ir_TM75)
plot(results_masked)

# mask and plot Jan1 temp predictions in IE
krg_results_df <- raster::as.data.frame(results_masked, xy = T)
krg_results_df <- krg_results_df[which(!is.na(krg_results_df$var1.pred)), ]
ggplot(data = krg_results_df, 
       aes(x = x, y = y, color = var1.pred)) +
  geom_point()


IE_krig_maps <- raster::brick(krg_predict)
IE_krig_maps <- raster::mask(IE_krig_maps, ir_TM75) # mask with coast shapefile
names(IE_krig_maps) <- c("Jan1  daily temp prediction", "variance")
plot(IE_krig_maps)

## !!---- RESULT ----!!
# At this point, krg_results_df is a useable result of interpolated temps for 
# Jan 1, 1995
## ----------------------- end interpolation ---------------------------------






##################### testing with smaller data ###############
# # testing how plots look
# # plot using lat/long
# small <- eobs_data[which(eobs_data$longitude < -8), ]
# ggplot(data = small, aes(x = longitude, y = latitude, color = elevation)) + 
#   geom_point()
# 
# # drop spatial class
# osi_df <- as.data.frame(osi)
# names(osi_df)[13:14] <- c("eastings", "northings")
# 
# ggplot(data = osi_df, aes(x = eastings, y = northings, 
#                                          color = tg)) + 
#   geom_point()
# 
# ## do the same thing for the lat/long version
# # drop spatial class
# spat_small_df <- as.data.frame(spat_small)
# 
# ggplot(data = spat_small_df, aes(x = longitude, y = latitude, 
#                                  color = tg)) + 
#   geom_point()
