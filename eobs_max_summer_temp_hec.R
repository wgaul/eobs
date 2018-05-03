#####################################
## Interpolate E-OBS maximum summer temp (98th quantile) to Irish hectad scale
## 
## inputs: eobs.RData - data prepared by Jon Yearsley and shared in the 
##            Grassland dropbox folder, summer 2017
## outputs: krg_summer_tx_predict - a spatialGridDataFrame holding 
##          interpolated values of 98th quantile of maximum summer temperature 
##          over the years 1995-2016 at the hectad scale.  
##          This includes many ocean blocks where 
##          predictions are likely poor.  This should be masked before being
##          used for analysis.  This is the object to use for most things.
##          
##          krg_summer_tx_rast - a rasterized version of krg_summer_tx_predict  
##          
##          krig_summer_tx_map - a rasterBrick version of krg_summer_tx_predict,
##          useful for easy plotting. krg_summer_tx_predict is probably the 
##          object to use for analysis.
## 
##          year_summer_tx - a list containing an element for each year 
##          1995-2016 and each element is a list of length 3 
##          holding the krig predictions as a spatial grid, as a raster, and
##          as a rasterBrick (for easy plotting)
##          
## NOTES: This script first finds the 98th quantile of maximum daily summer
## temperatures for each year, then takes the mean over all years to give 
## mean summer max temp.  It then interpolates that to the hectad scale.
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

## calculate 98th quantile of maximum daily temp for each year -----------------
# tx is maximum daily temperature (deg C). 
eobs_data$location <- paste(eobs_data$latitude, eobs_data$longitude, sep = "_")
eobs_data$year <- lubridate::year(eobs_data$date)
eobs_data$month <- lubridate::month(eobs_data$date)
tx_summary <-  select(eobs_data, longitude, latitude, 
                             location, month, year, tx) %>%
  filter(month < 10 & month > 3) %>%
  group_by(location, year) %>%
  mutate(annual_tx = stats::quantile(tx, probs = 0.98)) %>%
  select(-tx, -month) %>%
  distinct()

# calculate mean max summer temp over all years
tx_summary <- ungroup(tx_summary) %>%
  group_by(location) %>%
  mutate(mean_tx = mean(annual_tx, na.rm = T))

mean_tx <- select(tx_summary, 
                  -c(year, annual_tx)) %>%
  distinct()

## end calculating annual min winter temp -------------------------------------

## make both yearly and average yearly precip into spatial dfs -----------------
## yearly 
# make yearly summaries into spatial data frames for interpolation
spat_tx <- tx_summary[, 3:ncol(tx_summary)]
coordinates(spat_tx) <- as.matrix(tx_summary[, 1:2])
proj4string(spat_tx) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north -----
spat_tx <- spTransform(spat_tx, CRS("+init=epsg:29903"))
dimnames(spat_tx@coords)[[2]][which(
  dimnames(spat_tx@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_tx@coords)[[2]][which(
  dimnames(spat_tx@coords)[[2]] == "latitude")] <- "northings"

## mean annual
# make multi-year temp summaries into spatial data frames for interpolation
spat_mean_tx <- mean_tx[, 3:ncol(mean_tx)]
coordinates(spat_mean_tx) <- as.matrix(mean_tx[, 1:2])
proj4string(spat_mean_tx) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_mean_tx <- spTransform(spat_mean_tx, 
                            CRS("+init=epsg:29903"))
dimnames(spat_mean_tx@coords)[[2]][which(
  dimnames(spat_mean_tx@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_mean_tx@coords)[[2]][which(
  dimnames(spat_mean_tx@coords)[[2]] == "latitude")] <- "northings"

## ----------------------- interpolate by krigging ----------------------------
# Interpolate data in spat_tx and spat_mean_tx to the raster 
# grid irish_hec_raster
irish_spat_grid <- as(irish_hec_raster, 'SpatialGrid')

## mean min (98th quantile) winter temp over all years
# create empirical variogram for mean min winter temp
gs_mean_tx <- gstat(id = "mean_tx", 
                    formula = mean_tx ~ 1, 
                    data = spat_mean_tx) # make gstat object
v_mean_tx <- variogram(gs_mean_tx)
#vgm(as.character(vgm()[, 1])) # show vgm options
f_var_mean_tx <- fit.variogram(v_mean_tx, 
                               vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", 
                                     "Ste"))) 
f_var_mean_tx
plot(variogramLine(f_var_mean_tx, max(v_mean_tx$dist)), type = 'l')
points(v_mean_tx[, 2:3], pch = 20, col = 'red')

# use variogram in kriging interpolation
krg_mean_tx <- gstat(formula = mean_tx~1, 
                     data = spat_mean_tx, model = f_var_mean_tx)

krg_mean_tx_predict <- predict(krg_mean_tx, irish_spat_grid)
names(krg_mean_tx_predict) <- c("mean_tx", "variance")

# make rasters
krg_mean_tx_rast <- raster::raster(krg_mean_tx_predict)

# masking with shapefile, but I think the shapefile cuts off some hectads that
# do have data (some sp. datasets have 1014 hectads)
plot(raster::mask(krg_mean_tx_rast, ir_TM75))

krig_mean_tx_map <- raster::brick(krg_mean_tx_predict)
krig_mean_tx_map <- raster::mask(krig_mean_tx_map, ir_TM75) 
names(krig_mean_tx_map) <- c("mean_tx", "variance")
plot(krig_mean_tx_map, main = "Average 98th quantile of summer tx")

## make list of interpolated 98th quantile of winter tn for each year
interp <- function(year, spat_data, var_name) {
  # create empirical variogram for winter tn
  spat_data <- spat_data[spat_data$year == year, ]
  gs <- gstat(id = "summer_tx", formula = annual_tx~1, 
                 data = spat_data) # make gstat object
  v <- variogram(gs)
  f_var <- fit.variogram(v, vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", "Ste"))) 
  # plot(variogramLine(f_var, max(v$dist)), type = 'l', 
  #      ylim = c(0, max(v[, 3])))
  # points(v[, 2:3], pch = 20, col = 'red')
  
  # use variogram in kriging interpolation
  krg <- gstat(formula = annual_tx~1, 
               data = spat_data, model = f_var)
  
  krg_predict <- predict(krg, irish_spat_grid)
  names(krg_predict) <- c("summer_tx", "variance")
  
  # make rasters
  krg_rast <- raster::raster(krg_predict)
  krig_map <- raster::brick(krg_predict)
  krig_map <- raster::mask(krig_map, ir_TM75) 
  names(krig_map) <- c("summer_tx", "variance")
  
  krg_results <- list(hec_predict = krg_predict, 
                      hec_rast = krg_rast, 
                      hec_map = krig_map)
  krg_results
}

year_summer_tx <- lapply(unique(spat_tx$year), 
                         spat_data = spat_tx, 
                         FUN = interp)
names(year_summer_tx) <- unique(spat_tx$year)

# plot just to make sure it looks ok
pdf("./yearly_summer_tx_plots.pdf")
for (i in 1:length(year_summer_tx)) {
  print(plot(year_summer_tx[[i]]$hec_map, main = names(year_summer_tx)[i]))
}
dev.off()

### save outputs ---------------------------------------------------------------
save(krg_mean_tx_predict, krg_mean_tx_rast, krig_mean_tx_map, year_summer_tx, 
     file = "summer_tx_hectad.RData")