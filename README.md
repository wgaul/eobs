# eobs
Manipulating and Interpolating EObs climate data
author: Willson Gaul
created: October 2017
last modified: 8 Nov 2017

## Results of interpolation scripts:
The interpolation scripts produce objects that I have saved in .RData files.  Here are brief descriptions of those files:

- **annual_precip_hectad.RData** holds spatial dataframes and rasters produced by the "eobs_annual_precipitation_hectad.R" script.  These objects have mean annual precipitation over the years 1995-2016 (excluding 2010-2012) and also annual precipitation during each year, both interpolated from EObs data to the hectad scale.  
- **elevation_hec.RData** The spatial dataframe produced by the eobs_elevatoin_hectad.R script, which interpolates elevation from the EObs data to the hectad scale.  Note that the maximum altitudes produced by this krigging interpolation are much lower than the alts of actual peaks, probably because EObs data is recorded at a coarse scale.  I should probably re-do this using the original USGS data that the EObs elevations are based on, because I think the USGS data is at a finer scale and would more accurately represent peaks, especially once I get around to interpolating to the 1km grid scale.

