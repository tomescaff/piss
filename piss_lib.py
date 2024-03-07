#!/home/alvaro/.conda/envs/conda/bin/python
# -*- coding: utf-8 -*-
#
# ~/projects/piss/src/piss_lib.py
#
# info  : common functions for piss project.
# usage : ---
# author: @alvaroggc


# standard libraries
import locale
import warnings
from functools import partial
from glob import glob


# 3rd party packages
import pyproj
import numpy as np
import xarray as xr
import multiprocess as mp
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# function aliases
print = partial(print, flush=True)

# figures configuration
plt.rcParams['figure.dpi'     ] = 300
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8

# general configuration
warnings.filterwarnings('ignore')               # supress deprecation warnings
locale.setlocale(locale.LC_ALL, 'es_CL.UTF-8')  # apply spanish locale settings

# local source
from piss_params import *


############################
##### GLOBAL FUNCTIONS #####
############################


def log_message(msg):
    '''
    info: print log header message as "\n:: {string}: ---" and return indent parameter.
    parameters:
        msg : str -> header message to log.
    returns:
        indent : str -> string with spaces of length (len(string) + 5).
    '''

    # logging message
    print(f'\n:: {msg}: ---')

    # indent parameter
    indent = " " * (len(msg) + 5)

    # output
    return indent


def load_simulation(dirin, simid, vargroup='h1', keys=None):
    '''
    info:
    parameters:
    returns:
    '''

    # list of netcdf files to load
    if vargroup == 'h0':

        # main location
        fin = sorted(glob(f'{dirin}/cesm/archives/*.{simid.upper()}.*/atm/hist/*.cam.h0.*.nc'))

    elif vargroup == 'h1':

        # list of netcdf files to load
        fin = sorted(glob(f'{dirin}/data/piss/h1/yData/{simid.lower()}/*.nc'))

    # logging message
    indent = log_message('loading files')
    for f in fin: print(f"{indent}{f.split('/')[-1]}")

    # load all files inside dataset container
    ds = xr.open_mfdataset(fin)#, decode_times=False)

    # select variables
    if keys: ds = ds[keys]

    # round spatial coordinates
    if 'lat' in ds.coords.keys(): ds['lat'] = ds['lat'].where(False, np.round(ds['lat'], 5))
    if 'lon' in ds.coords.keys(): ds['lon'] = ds['lon'].where(False, np.round(ds['lon'], 5))

    # make longitude cyclic
    dslast        = ds.sel({'lon': 0}).copy()
    dslast['lon'] = 360
    ds            = xr.concat((ds, dslast), 'lon')

    # output
    return ds


def xarray_time_iterator(da):
    '''
    info:
    parameters:
    returns:
    '''

    # create iterator
    for t in da['time']:

        # generate iterator item
        yield da.sel({'time': t})


def convert_to_lambert(da, res=90):
    '''
    info:
    parameters:
    returns:
    '''

    def regrid_timestep(dak):
        '''
        info:
        parameters:
        returns:
        '''

        # project data
        dak_xy = griddata(points=(X.flatten(), Y.flatten()),
                          values=dak.data.flatten(),
                          xi=(XN, YN),
                          method='cubic')

        # output
        return dak_xy

    # create lat-lon grids
    LON, LAT = np.meshgrid(da['lon'].data, da['lat'].data)

    # create transformer
    latlon_to_lambert = pyproj.Transformer.from_crs(latlon_crs, laea_crs)

    # transform curvilinear grid (lat-lon) to cartesian grid (x-y)
    X, Y = latlon_to_lambert.transform(LAT, LON)

    # define new regular x-y grid
    lim    = np.floor(np.abs(X.max() / 100)) * 100
    xn     = np.arange(-lim, lim+1, res)
    yn     = np.arange(-lim, lim+1, res)
    XN, YN = np.meshgrid(xn, yn, indexing='xy')

    # create container of projected data
    newarr = np.zeros((len(da['time']), len(xn), len(yn)))

    # create thread pool
    pool = mp.Pool(processes=25)

    # create arguments iterator
    args = xarray_time_iterator(da)

    # compute processing
    results = pool.imap(regrid_timestep, args)

    # close thread pool
    pool.close()
    pool.join()

    #
    # process each timestep
    for k, dak_xy in enumerate(results):

        # add results to container
        newarr[k, ::] = dak_xy

    # make new container for regular x-y data
    nda = xr.DataArray(data=newarr,
                       dims=['time', 'y', 'x'],
                       coords={'time': da['time'], 'x': xn, 'y': yn},
                       attrs=da.attrs)

    # assign x-y grid coordinates to new data
    nda = nda.assign_coords({'X': (('y', 'x'), XN), 'Y': (('y', 'x'), YN)})

    # output
    return nda


def distance_between_points(p1, p2):
    '''
    info:
    parameters:
    returns:
    '''

    # separate coordinates
    x1 = p1[1]
    x2 = p2[1]
    y1 = p1[0]
    y2 = p2[0]

    # convert distance to km
    dist = np.sqrt((x1-x2)**2 + (y1-y2)**2)

    # output (minimum distance)
    return dist


def contour_search(da, point):
    '''
    info:
    parameters:
    returns:
    '''

    # separate coordinates of field
    lat   = da['lat']
    lon   = da['lon']
    dates = da['time']

    # separate coordinates of point
    latp  = point[0]
    lonp  = point[1]
    datep = cf.dt(point[-1], calendar='noleap')

    # find indexes of point inside coordinates
    ilat, = np.where(np.abs((lat - latp)) == np.min(np.abs((lat - latp))))[0]
    jlon, = np.where(np.abs((lon - lonp)) == np.min(np.abs((lon - lonp))))[0]
    kt  , = np.where(np.abs((dates - datep)) == np.min(np.abs((dates - datep))))[0]

    # choose date
    # TODO: dak = da.isel({'time': it}).squeeze()

    # process each contour



