#!/home/cr2/agomez/.conda/envs/conda/bin/python
# -*- coding: utf-8 -*-
#
# ~/cesm/lib/thesis_parameters
#
# info  : common parameters used for processing CESM simulations.
# author: @alvaroggc


# 3rd party packages
import cf
import pyproj
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


#############################
##### GLOBAL PARAMETERS #####
#############################


# range of dates to process
date_ini = '0001-01-01 00:00:00'
date_end = '0001-12-31 23:59:59'

# projection definitions
latlon_crs = 'epsg:4326'
laea_crs   = '+proj=laea +lon_0=-73.53 +lat_0=-90 +ellps=WGS84 +x_0=0 +y_0=0 +units=km'

# map definition
proj  = ccrs.NearsidePerspective(central_longitude=(360-73.53), central_latitude=-90)
trans = ccrs.PlateCarree()

# datetime normalization parameters
ref_datetime = cf.dt(2000, 1, 1, calendar='gregorian')
norm         = np.timedelta64(1, 'D')
