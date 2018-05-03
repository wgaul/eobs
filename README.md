# eobs

Manipulating and Interpolating EObs climate data

created: October 2017

last modified: 3 May 2018

## Results of interpolation scripts:
The interpolation R scripts produce objects that get saved in .RData files.  Here are brief descriptions of the outputs created by each R script:

- **eobs_annual_precipitation_hectad.R** holds spatial dataframes and rasters produced by the "eobs_annual_precipitation_hectad.R" script.  These objects have mean annual precipitation over the years 1995-2016 (excluding 2010-2012) and also annual precipitation during each year, both interpolated from EObs data to the hectad scale.  

- **mean_pp_hectad.RData** 



Data used here are daily gridded data (0.25 deg resolution) for Ireland (1995-2016). The data are created as part of the E-OBS gridded data product from the European Climate Assessment and Dataset (ECA&D) EU project http://www.ecad.eu/download/ensembles/ensembles.php
 
Data are downloaded from http://www.ecad.eu/download/ensembles/downloadchunks.php

Variables in this data set are:

- tg = daily mean temperature (units = deg Celsius)
- tn = daily minimum temperature (units = deg Celsius)
- tx = daily maximum temperature (units = deg Celsius)
- rr = daily precipitation sum (units = mm)
- pp = daily averaged sea level pressure  (units = hecto Pascals)
- elevation = mean elevation (m)

Each variable has a standard error

Dataset details: 

http://www.ecad.eu/download/ensembles/Haylock_et_al_2008.pdf

http://onlinelibrary.wiley.com/doi/10.1029/2010JD015468/full