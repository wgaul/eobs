# eobs

Manipulating and Interpolating E-OBS climate data

These scripts interpolate E-OBS climate data to the Irish 10-km squares (hectads) on the Irish OSI grid system (epsg 29903).  

## About E-OBS data
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