#####################################
## Interpolate E-OBS minimum winter temp (2nd quantile) to Irish hectad scale
## 
## inputs:  * eobs.RData - data prepared by Jon Yearsley and shared in the 
##            Grassland dropbox folder, summer 2017
##          * ireland_coastline shapefile (I got this from Jon Yearsley but I
##                think it might be a naturalEarth file originally)
##          * IE_10km_hecs.shp shapefile of Irish hectads created in QGIS by wg
## outputs: krg_min_winter_tmp_predict - a spatialGridDataFrame holding 
##          interpolated values of 2nd quantile of minimum winter temperature 
##          over the years 1995-2016 at the hectad scale.  
##          This includes many ocean blocks where 
##          predictions are likely poor.  This should be masked before being
##          used for analysis.  This is the object to use for most things.
##          
##          krg_min_winter_tmp_rast - a rasterized version of 
##          krg_min_winter_tmp_predict  
##          
##          krig_min_winter_tmp_map - a rasterBrick version of 
##          krg_min_winter_tmp_predict, useful for easy plotting.  
##          krg_min_winter_tmp_predict is probably the object to use for 
##          analysis.
## 
##          year_min_winter_tmp - a list containing an element for each year 
##          1995-2016 and each element is a list of length 3 
##          holding the krig predictions as a spatial grid, as a raster, and
##          as a rasterBrick (for easy plotting)
##          
## NOTES: This script first finds the 2nd quantile of minimum daily winter
## temperatures for each year, then takes the mean over all years to give 
## mean min winter temp.  It then interpolates that to the hectad scale.
##
## TODO:  - Perhaps try the pseudo-matlab version Jon mentioned?
##          alternatviely, try approxfun() or the R equivalent of Matlab 
##          "interp" "interp2" for 2d.
## 
## author: Willson Gaul
## created: 13 Dec 2017
## last modified: 13 Dec 2017
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

## calculate 2nd quantile of minimum daily temp for each year -----------------
# tn is minimum daily temperature (deg C). 
eobs_data$location <- paste(eobs_data$latitude, eobs_data$longitude, sep = "_")
eobs_data$year <- lubridate::year(eobs_data$date)
eobs_data$month <- lubridate::month(eobs_data$date)
min_tmp_summary <-  select(eobs_data, longitude, latitude, 
                           location, month, year, tn) %>%
  filter(month >= 9 | month <= 4) %>%
  group_by(location, year) %>%
  mutate(annual_min_winter_tmp = stats::quantile(tn, probs = 0.02)) %>%
  select(-tn, -month) %>%
  distinct()

# calculate mean minimum temp over all years
min_tmp_summary <- ungroup(min_tmp_summary) %>%
  group_by(location) %>%
  mutate(mean_annual_min_winter_tmp = mean(annual_min_winter_tmp, na.rm = T))

mean_annual_min_winter_tmp <- select(min_tmp_summary, 
                                     -c(year, annual_min_winter_tmp)) %>%
  distinct()

## end calculating annual min winter temp -------------------------------------

## make both yearly and average yearly precip into spatial dfs -----------------
## yearly 
# make yearly summaries into spatial data frames for interpolation
spat_min_temp <- min_tmp_summary[, 3:ncol(min_tmp_summary)]
coordinates(spat_min_temp) <- as.matrix(min_tmp_summary[, 1:2])
# set proj4string based on WGS84 here: http://spatialreference.org/ref/epsg/4326/
proj4string(spat_min_temp) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north -----
spat_min_temp <- spTransform(spat_min_temp, CRS("+init=epsg:29903"))
dimnames(spat_min_temp@coords)[[2]][which(
  dimnames(spat_min_temp@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_min_temp@coords)[[2]][which(
  dimnames(spat_min_temp@coords)[[2]] == "latitude")] <- "northings"

## mean annual
# make multi-year temp summaries into spatial data frames for interpolation
spat_mean_tn <- mean_annual_min_winter_tmp[, 3:ncol(
  mean_annual_min_winter_tmp)]
coordinates(spat_mean_tn) <- as.matrix(
  mean_annual_min_winter_tmp[, 1:2])
# set proj4string based on WGS84 here: http://spatialreference.org/ref/epsg/4326/
proj4string(spat_mean_tn) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_mean_tn <- spTransform(spat_mean_tn, 
                                     CRS("+init=epsg:29903"))
dimnames(spat_mean_tn@coords)[[2]][which(
  dimnames(spat_mean_tn@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_mean_tn@coords)[[2]][which(
  dimnames(spat_mean_tn@coords)[[2]] == "latitude")] <- "northings"
# end coordinate conversion -----

## ----------------------- interpolate by krigging ----------------------------
# Interpolate data in spat_min_temp and spat_allyear_min_tempto the raster 
# grid irish_hec_raster
irish_spat_grid <- as(irish_hec_raster, 'SpatialGrid')

## mean min (2nd quantile) winter temp over all years
# create empirical variogram for mean min winter temp
gs_mean_tn <- gstat(id = "mean_min_winter_temp", 
                    formula = mean_annual_min_winter_tmp ~ 1, 
                    data = spat_mean_tn) # make gstat object
v_mean_tn <- variogram(gs_mean_tn)
#vgm(as.character(vgm()[, 1])) # show vgm options
f_var_mean_tn <- fit.variogram(v_mean_tn, 
                               vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", 
                                     "Ste"))) 
f_var_mean_tn
plot(variogramLine(f_var_mean_tn, max(v_mean_tn$dist)), type = 'l')
points(v_mean_tn[, 2:3], pch = 20, col = 'red')

# use variogram in kriging interpolation
krg_mean_tn <- gstat(formula = mean_annual_min_winter_tmp~1, 
                     data = spat_mean_tn, model = f_var_mean_tn)

krg_mean_tn_predict <- predict(krg_mean_tn, irish_spat_grid)
names(krg_mean_tn_predict) <- c("mean_winter_tn", "variance")

# make rasters
krg_mean_tn_rast <- raster::raster(krg_mean_tn_predict)

# masking with shapefile, but I think the shapefile cuts off some hectads that
# do have data (some sp. datasets have 1014 hectads)
plot(raster::mask(krg_mean_tn_rast, ir_TM75))

krig_mean_tn_map <- raster::brick(krg_mean_tn_predict)
krig_mean_tn_map <- raster::mask(krig_mean_tn_map, ir_TM75) 
names(krig_mean_tn_map) <- c("mean_tn", "variance")
plot(krig_mean_tn_map, main = "Average 2nd quantile of winter tn")


## make list of interpolated 2nd quantile of winter tn for each year
interp <- function(year, spat_data, var_name) {
  # create empirical variogram for winter tn
  spat_data <- spat_data[spat_data$year == year, ]
  gs_tn <- gstat(id = "winter_tn", formula = annual_min_winter_tmp~1, 
                 data = spat_data) # make gstat object
  v <- variogram(gs_tn)
  f_var <- fit.variogram(v, vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", "Ste"))) 
  # plot(variogramLine(f_var, max(v$dist)), type = 'l', 
  #      ylim = c(0, max(v[, 3])))
  # points(v[, 2:3], pch = 20, col = 'red')
  
  # use variogram in kriging interpolation
  krg <- gstat(formula = annual_min_winter_tmp~1, 
               data = spat_data, model = f_var)
  
  krg_predict <- predict(krg, irish_spat_grid)
  names(krg_predict) <- c("winter_tn", "variance")
  
  # make rasters
  krg_rast <- raster::raster(krg_predict)
  krig_map <- raster::brick(krg_predict)
  krig_map <- raster::mask(krig_map, ir_TM75) 
  names(krig_map) <- c("winter_tn", "variance")
  
  krg_results <- list(hec_predict = krg_predict, 
                      hec_rast = krg_rast, 
                      hec_map = krig_map)
  krg_results
}

year_winter_tn <- lapply(unique(spat_min_temp$year), 
                         spat_data = spat_min_temp, 
                         FUN = interp)
names(year_winter_tn) <- unique(spat_min_temp$year)

# plot just to make sure it looks ok
pdf("./yearly_winter_tn_plots.pdf")
for (i in 1:length(year_winter_tn)) {
  print(plot(year_winter_tn[[i]]$hec_map, main = names(year_winter_tn)[i]))
}
dev.off()

### save outputs ---------------------------------------------------------------
save(krg_mean_tn_predict, krg_mean_tn_rast, krig_mean_tn_map, year_winter_tn, 
     file = "winter_tn_hectad.RData")
