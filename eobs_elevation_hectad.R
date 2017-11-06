###################################
## Interpolate E-OBS altitude data (1995-2016) to Irish hectad scale
## 
## inputs: eobs.RData - data prepared by Jon Yearsley and shared in the 
##            Grassland dropbox folder, summer 2017
## outputs: krg_elev_predict - a spatial data frame of elevation at the 
##              location and scale of Irish hectads.  This object has not been
##              masked or trimmed to fit Ireland, so it has lots of cells over
##              the ocean for which predictions are probably not good.  
##              However, I don't want to mask this or filter it yet because 
##              the IE coastline shapefile excludes some hectads that have 
##              land and data.  So this will need to be subset in scripts that
##              use it for analysis.
##              
## NOTES: This drastically under estimates elevations of peaks, probably because 
## those areas don't stay high over the large area that eobs elevation data
## is recorded at.  I think the eobs elevation data came from some USGS source,
## which might have finer resolution, so I should probably see if I can 
## interpolate elevation to hectads straight from the original source, rather 
## than from eobs, which has probably been through one interpolation already.
## 
## author: Willson Gaul
## created: 25 Sep 2017
## last modified: 6 Nov 2017
#################################### 

setwd("~/Documents/Data_Analysis/UCD/predictor_variables/eobs")
library(wgutil)
library(rgdal)
library(gdalUtils)
library(ggplot2)
library(fields)
library(gstat)
library(raster)
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

## ------------------ convert eobs lat/long to OSI east/north ----------------
# make eobs_data a spatial data frame
spat_eobs <- eobs_data[, 3:14]
coordinates(spat_eobs) <- as.matrix(eobs_data[, 1:2])
# set proj4string based on WGS84 here: http://spatialreference.org/ref/epsg/4326/
proj4string(spat_eobs) <- CRS(
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# convert lat/long to OSI east/north
spat_eobs <- spTransform(spat_eobs, CRS("+init=epsg:29903"))
dimnames(spat_eobs@coords)[[2]][which(
  dimnames(spat_eobs@coords)[[2]] == "longitude")] <- "eastings"
dimnames(spat_eobs@coords)[[2]][which(
  dimnames(spat_eobs@coords)[[2]] == "latitude")] <- "northings"

## -------------------------- end coordinate conversion ------------------------

## ------------------------- interpolate kriggin ------------------------------
# I want to interpolate data in spat_eobs to the raster grid irish_hec_raster
# following this: http://rspatial.org/analysis/rst/4-interpolation.html
irish_spat_grid <- as(irish_hec_raster, 'SpatialGrid')

# create empirical variogram for elevation
# Elevation doesn't change over time, so can use values for just 1 day
elev <- spat_eobs[spat_eobs$date == as.character(median(spat_eobs$date)), ] 
gs_eobs <- gstat(id = "elevation", formula = elevation~1, 
                 data = elev) # make gstat object
v <- variogram(gs_eobs, width = 1000)

#vgm(as.character(vgm()[, 1])) # show vgm options
f_var <- fit.variogram(v, vgm(c("Exp", "Sph", "Gau", "Exc", "Mat", 
                                "Ste"))) 

# alternatviely, try approxfun() or the R equivalent of Matlab "interp" "interp2" for 2d.

f_var
plot(variogramLine(f_var, max(v$dist)), type = 'l', 
     ylim = c(0, max(v[, 3])))
points(v[, 2:3], pch = 20, col = 'red')

# use variogram in kriging interpolation
krg <- gstat(formula = elevation~1, data = elev, model = f_var)

krg_elev_predict <- predict(krg, irish_spat_grid)

## !!! ---- RESULT -------- !!!
## krg_elev_predict has the results I want
krg_elev_rast <- raster::raster(krg_elev_predict)

# masking with shapefile, but I think the shapefile cuts off some hectads that
# do have data (some sp. datasets have 1014 hectads)
plot(raster::mask(krg_elev_rast, ir_TM75))

IE_krig_maps <- raster::brick(krg_elev_predict)
IE_krig_maps <- raster::mask(IE_krig_maps, ir_TM75) # mask with coast shapefile
names(IE_krig_maps) <- c("Elevation", "variance")
plot(IE_krig_maps)
## ----------------------- end interpolation ---------------------------------




## ---------------------- linear interpolation? -----------------------------
# I think these methods will first interpolate a surface using original points,
# then subsequently rasterize that to hectad scale (so averaging to hectad 
# blocks happens after interpolation).  Don't know how the kriging method does
# this.
if(do_linear == T){
  
  # Elevation doesn't change over time, so can use values for just 1 day
  elev <- spat_eobs[spat_eobs$date == as.character(median(spat_eobs$date)), ] 
  
  # set mean of all observations as baseline to compare methods to
  RMSE <- function(observed, predicted) {
    sqrt(mean((predicted - observed)^2, na.rm = T))
  }
  # get RMSE for null model
  null <- RMSE(elev@data$elevation, 
               predicted = mean(elev@data$elevation))
  
  # control points
  cp <- mask(irish_hec_raster, ir_TM75)
  cp <- rasterToPoints(cp)
  
  # distance matrix
  d <- pointDistance(cp[, 1:2], elev@coords, lonlat = F)
  
  # find 4 nearest neighbours to each point
  nn <- 6
  nearest <- t(apply(d, 1, function(x) order(x)[1:nn]))
  
  # check if this looks right
  plot(ir_TM75)
  points(cp[601, 1:2, drop=F], col = 'blue', cex = 1)
  points(elev[nearest[601, ], ], col = 'red')
  
  # make pairs
  pairs <- cbind(rep(1:nrow(nearest), nn), as.vector(nearest))
  
  values <- elev$elevation[pairs[, 2]]
  lin_interp <- tapply(values, pairs[, 1], mean)
  
  # make raster to hold interpolated results
  nn_raster <- mask(irish_hec_raster, ir_TM75)
  # assign interpolated values to raster cells
  nn_raster[!is.na(nn_raster)] <- lin_interp 
  nn_points <- rasterToPoints(nn_raster)
  plot(nn_raster)
  
  # rmse using all data
  RMSE(observed = elev$elevation, 
       predicted = raster::extract(nn_raster, elev, method = "simple"))
  
  # cross-validate result
  library(dismo)
  kf <- kfold(nrow(elev))
  rast_dup <- mask(irish_hec_raster, ir_TM75)
  nn_rmse <- rep(NA, 5)
  for (k in 1:5) {
    test <- elev[kf == k, ]
    train <- elev[kf != k, ]
    d <- pointDistance(cp[, 1:2], train, lonlat = F)
    ngb <- t(apply(d, 1, function(x) order(x)[1:nn]))
    pairs <- cbind(rep(1:nrow(ngb), nn), as.vector(ngb))
    values <- elev$elevation[pairs[, 2]]
    pn <- tapply(values, pairs[, 1], mean)
    rast_dup[!is.na(rast_dup)] <- pn
    p <- extract(rast_dup, test)
    nn_rmse[k] <- RMSE(test$elevation, p)
  }
  
  # compare rmse from different methods
  mean(nn_rmse)
  null
  1 - (mean(nn_rmse) / null)
  
  k_rmse <- RMSE(observed = elev$elevation, 
                 predicted = raster::extract(krg_elev_rast, 
                                             elev, 
                                             method = "simple"))
  null
  mean(nn_rmse)
  k_rmse
  # krigging appears much better (if I'm calculating this correctly)
  
  # kriggin gives a wider range of elevation values
  # both methods appear to give qualitatively similar results, but krigging gives
  # higher mountains
  # Given rmse results and the maps, I prefer the krigging method to this 
  # nearest neighbor linear method
  plot(mask(krg_elev_rast, ir_TM75), main = "krig")
  plot(nn_raster, main = "linear")
  plot(mask(krg_elev_rast, ir_TM75) - nn_raster)
  
} # end if(do_linear)