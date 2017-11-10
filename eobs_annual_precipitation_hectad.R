###################################
## Interpolate E-OBS total annual precipitation (1995-2016) to Irish hectad scale
## 
## inputs: eobs.RData - data prepared by Jon Yearsley and shared in the 
##            Grassland dropbox folder, summer 2017
## outputs: krg_mean_rr_predict - a spatialGridDataFrame holding interpolated
##          values of mean annual precipitation over the years 1995-2016 
##          (excluding 2010-2012) at the hectad scale.  
##          This includes many ocean blocks where 
##          predictions are likely poor.  This should be masked before being
##          used for analysis.  This is the object to use for most things.
##          
##          krg_mean_rr_rast - a rasterized version of krg_mean_rr_predict.  
##          
##          krig_mean_rr_map - a rasterBrick version of krg_mean_rr_predict, 
##          useful for easy plotting.  krg_mean_rr_predict is probably the 
##          object to use for analysis.
## 
##          year_tot_rr - a list containing an element for each year 1995-2016
##          (excluding 2010-2012), and each element is a list of length 3 
##          holding the krig predictions as a spatial grid, as a raster, and
##          as a rasterBrick (for easy plotting)
##          
##              
## NOTES: This script first summarizes daily precipitation to produce annual
## precipitation amounts, then interpolates that to the hectad scale. 
##
## TODO:  - linear interpolation
##        - compare krig and linear results
##        - Perhaps try the pseudo-matlab version Jon mentioned?
# alternatviely, try approxfun() or the R equivalent of Matlab "interp" "interp2" for 2d.
## 
## author: Willson Gaul
## created: 25 Sep 2017
## last modified: 8 Nov 2017
#################################### 

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

do_linear <- F # do linear interpolation also (at bottom script)?

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

### calculate annual precipitation --------------------------------------------
# rr is daily precipitation in mm 
# There are some NA values in rr.  Diagnostic code at the end of this script
# looks at these.
# It seems all the NA values are from 2010-2012, and affect many hectads 
# (the majority of hectads in the Republic).  So, as an initial way of dealing
# with this, I propose just not using 2010-2012 in calculating average annual
# precipitation.

## calculate average annual precipitation from 1995-2015 (excluding 2010-2012)
eobs_data$location <- paste(eobs_data$latitude, eobs_data$longitude, sep = "_")
eobs_data$year <- lubridate::year(eobs_data$date)
precip_summary <-  select(eobs_data, longitude, latitude, 
                          location, year, rr) %>%
  group_by(location, year) %>%
  mutate(annual_rr_mm = sum(rr)) %>%
  select(-rr) %>%
  distinct()

# remove sums for the years that have NAs
drop_yrs <- unique(eobs_data$year[which(is.na(eobs_data$rr))])
precip_summary$annual_rr_mm[which(precip_summary$year %in% drop_yrs)] <- NA

## calculate mean annual precipitation
precip_summary <- ungroup(precip_summary) %>%
  group_by(location) %>%
  mutate(mean_annual_rr_mm = mean(annual_rr_mm, na.rm = T))

mean_annual_rr <- select(precip_summary, -c(year, annual_rr_mm)) %>%
  distinct()

## end calculating annual precipitation --------------------------------------

## make both yearly and average yearly precip into spatial dfs -----------------
## ------------------ convert eobs lat/long to OSI east/north ----------------
## yearly
# make precipitation summaries into spatial data frames for interpolation
spat_precip <- precip_summary[, 3:ncol(precip_summary)]
coordinates(spat_precip) <- as.matrix(precip_summary[, 1:2])
# set proj4string based on WGS84 here: http://spatialreference.org/ref/epsg/4326/
proj4string(spat_precip) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_precip <- spTransform(spat_precip, CRS("+init=epsg:29903"))
dimnames(spat_precip@coords)[[2]][which(
  dimnames(spat_precip@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_precip@coords)[[2]][which(
  dimnames(spat_precip@coords)[[2]] == "latitude")] <- "northings"

## mean annual
# make precipitation summaries into spatial data frames for interpolation
spat_mean_rr <- mean_annual_rr[, 3:ncol(mean_annual_rr)]
coordinates(spat_mean_rr) <- as.matrix(mean_annual_rr[, 1:2])
# set proj4string based on WGS84 here: http://spatialreference.org/ref/epsg/4326/
proj4string(spat_mean_rr) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_mean_rr <- spTransform(spat_mean_rr, CRS("+init=epsg:29903"))
dimnames(spat_mean_rr@coords)[[2]][which(
  dimnames(spat_mean_rr@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_mean_rr@coords)[[2]][which(
  dimnames(spat_mean_rr@coords)[[2]] == "latitude")] <- "northings"
## -------------------------- end coordinate conversion ------------------------

## ----------------------- interpolate by krigging ----------------------------
# I want to interpolate data in spat_eobs to the raster grid irish_hec_raster
# following this: http://rspatial.org/analysis/rst/4-interpolation.html
irish_spat_grid <- as(irish_hec_raster, 'SpatialGrid')

## mean annual precipitation
# create empirical variogram for mean annual precipitaion
gs_mean_rr <- gstat(id = "mean_yrl_rr", formula = mean_annual_rr_mm~1, 
                 data = spat_mean_rr) # make gstat object
v_mean_rr <- variogram(gs_mean_rr, width = 1000)
#vgm(as.character(vgm()[, 1])) # show vgm options
f_var_mean_rr <- fit.variogram(v_mean_rr, 
                               vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", 
                                "Ste"))) 
f_var_mean_rr
plot(variogramLine(f_var_mean_rr, max(v_mean_rr$dist)), type = 'l', 
     ylim = c(0, max(v_mean_rr[, 3])))
points(v_mean_rr[, 2:3], pch = 20, col = 'red')

# use variogram in kriging interpolation
krg_mean_rr <- gstat(formula = mean_annual_rr_mm~1, 
                     data = spat_mean_rr, model = f_var_mean_rr)

krg_mean_rr_predict <- predict(krg_mean_rr, irish_spat_grid)
names(krg_mean_rr_predict) <- c("mean_annual_rr", "variance")
## !!! ---- RESULT -------- !!!
## krg_mean_rr_predict has the results for mean annual precipitation from 
# 1995-2016 (excluding 2010-2012).  
# make rasters
krg_mean_rr_rast <- raster::raster(krg_mean_rr_predict)

# masking with shapefile, but I think the shapefile cuts off some hectads that
# do have data (some sp. datasets have 1014 hectads)
plot(raster::mask(krg_mean_rr_rast, ir_TM75))

krig_mean_rr_map <- raster::brick(krg_mean_rr_predict)
krig_mean_rr_map <- raster::mask(krig_mean_rr_map, ir_TM75) 
names(krig_mean_rr_map) <- c("mean_annual_rr", "variance")
plot(krig_mean_rr_map, main = "Average annual rr")

## make list of interpolated total rr for each year
interp <- function(year, spat_data, var_name) {
  # create empirical variogram for total precipitaion
  spat_data <- spat_data[spat_data$year == year, ]
  gs_rr <- gstat(id = "total_rr", formula = annual_rr_mm~1, 
                      data = spat_data) # make gstat object
  v <- variogram(gs_rr, width = 1000)
  #vgm() # show vgm options
  f_var <- fit.variogram(v, 
                         vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", 
                               "Ste"))) 
  # f_var
  # plot(variogramLine(f_var, max(v$dist)), type = 'l', 
  #      ylim = c(0, max(v[, 3])))
  # points(v[, 2:3], pch = 20, col = 'red')
  
  # use variogram in kriging interpolation
  krg <- gstat(formula = annual_rr_mm~1, 
               data = spat_data, model = f_var)
  
  krg_predict <- predict(krg, irish_spat_grid)
  names(krg_predict) <- c("total_rr", "variance")
  ## !!! ---- RESULT -------- !!!
  ## krg_predict has the results for total precipitation for this year 
  # make rasters
  krg_rast <- raster::raster(krg_predict)
  
  # masking with shapefile, but I think the shapefile cuts off some hectads that
  # do have data (some sp. datasets have 1014 hectads)
  #plot(raster::mask(krg_rast, ir_TM75))
  
  krig_map <- raster::brick(krg_predict)
  krig_map <- raster::mask(krig_map, ir_TM75) 
  names(krig_map) <- c("total_rr", "variance")
  #plot(krig_map, main = "Total annual rr")
  
  krg_results <- list(hec_predict = krg_predict, 
                      hec_rast = krg_rast, 
                      hec_map = krig_map)
  krg_results
}

keep_years <- unique(eobs_data$year)[which(
  unique(eobs_data$year) %nin% drop_yrs)]
names(keep_years) <- as.character(keep_years)
year_tot_rr <- lapply(keep_years, spat_data = spat_precip, 
                      FUN = interp)

# plot just to make sure it looks ok
pdf("./yearly_precipitation_plots.pdf")
for (i in 1:length(year_tot_rr)) {
  print(plot(year_tot_rr[[i]]$hec_map, main = names(year_tot_rr)[i]))
}
dev.off()

pdf("./yearly_precip_diff.pdf")
for (i in 2:(length(year_tot_rr))) {
  print(plot(year_tot_rr[[i]]$hec_map - year_tot_rr[[i-1]]$hec_map, 
             main = paste("Diff", names(year_tot_rr)[i])))
}
dev.off()




### Check NAs and errors ------------------------------------------------------
# these plots only work after creating spat_eobs based on 'nas' rather than
# based on eobs_data.  Don't need to view every time, diagnostic only
view_precip_nas <- F
if(view_precip_nas == T) {
  nas <- eobs_data[which(is.na(eobs_data$rr)), ]
  nas$location <- paste0(nas$latitude, nas$longitude)
  nas$year <- lubridate::year(nas$date)
  
  for (i in 1:length(unique(nas$location))) {
    print(unique(nas$location)[i])
    print(table(nas$year[nas$location == unique(nas$location)[i]]))
  }
  
  plot(ir_TM75)
  plot(spat_eobs[, ], add = T)
  plot(spat_eobs[which(spat_eobs$location == "51.875-8.125"), ], col = "red", 
       add = T)
  plot(spat_eobs[which(spat_eobs$location == "51.875-8.375"), ], col = "red", 
       add = T)
  plot(spat_eobs[which(spat_eobs$location == "53.875-9.375"), ], col = "red", 
       add = T)
}
## end diagnosing NAs

# check sums for site/year combinations with a for loop for a few locations \
# and years.
# This looks ok.
check_sums <- F
if(check_sums) {
  locs <- c("54.125_-7.375", "55.375_-5.625", "52.125_-7.875", "51.625_-9.625")
  check_df <- data.frame(matrix(ncol = 3, nrow = 100))
  colnames(check_df) <- c("location", "year", "sum_rr")
  r <- 1
  for (i in locs) {
    df <- eobs_data[eobs_data$location == i, ]
    df <- df[df$year %nin% c(2010, 2011, 2012), ]
    for (j in unique(df$year)) {
      check_df[r, "location"] <- i
      check_df[r, "year"] <- j
      check_df[r, "sum_rr"] <- sum(df[df$year == j, "rr"])
      r <- r + 1
    }
  }
}


### save outputs ---------------------------------------------------------------
save(krg_mean_rr_predict, krg_mean_rr_rast, krig_mean_rr_map, year_tot_rr, 
     file = "annual_precip_hectad.RData")
