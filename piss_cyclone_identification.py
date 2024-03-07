#!/home/cr2/agomez/.conda/envs/conda/bin/python
# -*- coding: utf-8 -*-
#
# ~/projects/piss/src/piss_cyclone_identification.py
#
# info  : identify cyclone centers in the souther hemisphere in piss simulations (cesm2), based on
#         Hanley & Caballero (2012).
# author: @alvaroggc


# standard libraries
import os
import re
import csv
import sys
import copy
import pickle
import datetime as dt

# 3rd party packages
import cf
import cartopy.crs as ccrs
import multiprocess as mp
from pyproj import Geod

# local source
from piss_lib import *


###########################
##### LOCAL FUNCTIONS #####
###########################


def convert_point_to_latlon(point):
    '''
    info:
    parameters:
    returns:
    '''

    # create transformer
    lambert_to_latlon = pyproj.Transformer.from_crs(laea_crs, latlon_crs)

    # separate coordinates and value
    x = point[1]
    y = point[0]

    # project to lat-lon
    lat, lon = lambert_to_latlon.transform(x, y)

    # transform lon range (from [-180, 180] to [0, 360])
    lon = ((lon + 360) % 360)

    # output
    return [lat, lon]


def convert_point_to_xy(point):
    '''
    info:
    parameters:
    returns:
    '''

    # create transformer
    latlon_to_lambert = pyproj.Transformer.from_crs(latlon_crs, laea_crs)

    # separate coordinates and value
    lon = point[1]
    lat = point[0]

    # transform point to cartesian grid (x-y)
    x, y = latlon_to_lambert.transform(lat, lon)

    # output
    return [y, x]


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

    # calculate distance
    dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

    # output (minimum distance)
    return dist


def east_or_not(p1, p2, D=500):
    '''
    info:
    parameters:
    returns:
    '''

    # convert to lat-lon coordinate system
    lat1, lon1 = convert_point_to_latlon(p1)
    lat2, lon2 = convert_point_to_latlon(p2)

    # longitude alternative
    lon1aux = (lon1 - 360) if (lon1 > 180) else lon1
    lon2aux = (lon2 - 360) if (lon2 > 180) else lon2

    # check if p2 is east of p1
    east = False if ((lon2 < lon1) and (lon2aux < lon1aux)) else True

    # if not moving east, at least check if p2 doesn't move much from p1
    if not east:

        # distance between points
        dist = distance_between_points(p1, p2)

        # check if distance is low
        if (dist < D):

            # update flag to indicate that p2 is a valid point
            east = True

    # output
    return east


def guess_point(p1, p2):
    '''
    info:
    parameters:
    returns:
    '''

    # convert to lat-lon coordinate system
    lat1, lon1 = convert_point_to_latlon(p1)
    lat2, lon2 = convert_point_to_latlon(p2)

    # if point to close to 0°lon, change range of longitude
    if (lon1 > 320) or (lon2 > 320):

        # shift coordinates
        lon1 = (lon1 - 360) if (lon1 > 180) else lon1
        lon2 = (lon2 - 360) if (lon2 > 180) else lon2

    # guess next point
    lat3 = lat2 + 0.75 * (lat2 - lat1)
    lon3 = lon2 + 0.75 * (lon2 - lon1)

    # convert to x-y coordinate
    p_guess = convert_point_to_xy([lat3, lon3])

    # output
    return p_guess


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


def get_variables(args: "list[str]"):
    '''
    info: retrieve variables from input arguments.
    parameters:
        args : list[str] -> list of arguments to parse.
    returns:
        simid : str  -> list of simulation(s) name(s) (default is 'lgm_100').
        plot  : bool -> plot maps with cyclones
        calc  : bool -> search for cyclones
    '''

    # get copy of arguments
    args = args[1:].copy()

    # formats
    simid_fmt = re.compile(r'^lgm_\d{3}?')      # simulation id
    D_fmt     = re.compile(r'^D\d')      # search distance

    # retrieve variables
    simid = [arg for arg in args if simid_fmt.match(str(arg))]   # sim. ID
    # D     = [arg for arg in args if D_fmt.match(    str(arg))]   # sim. ID
    calc  = True if 'calc' in args else False
    plot  = True if 'plot' in args else False

    # check arguments
    simid = 'lgm_100' if not simid else simid[0]
    # D     = 1000      if not D     else int(D[0][1:])

    # output
    return (simid, calc, plot)


###########################
##### LOCAL VARIABLES #####
###########################


# get variables from input arguments
simid, calc, plot = get_variables(sys.argv)


############################
##### LOCAL PARAMETERS #####
############################


# directories
homedir = os.path.expanduser('~')               # home directory
dirin   = f'/mnt/cirrus/results/friquelme'      # cesm simulations
dirout  = f'{homedir}/projects/piss/data'       # data output
dirimg  = f'{homedir}/projects/piss/img'        # output for figures

# output file with cyclones information
fout = f'cyclones_{simid}.pkl'

# variable to process
key = 'PSL'

# indexer to choose southern hemisphere (from -30° to the south)
shidx = {'lat': slice(-90, -5)}

# parameters needed in method
topo_max = 1000     # slp minimum points over this altitude are eliminated [m]
radius   = 1000     # gradient search radius [km]
grad_min = 10       # minimum slp gradient required for cyclone [hPa]

# figure parameters
size = (7.5, 6)


###################
##### RUNNING #####
###################


# first check if data should be loaded
if ('slp_xy' not in locals()):

    # load simulation datasets
    ds1 = load_simulation(dirin, simid, 'h1', ['PSL' ])
    ds0 = load_simulation(dirin, simid, 'h0', ['PHIS'])

    # separate variables
    slp  = ds1['PSL' ]
    topo = ds0['PHIS']

    # extract temporal range (only for slp)
    slp  = slp.sel( {'time': slice(date_ini, date_end)})

    # extract southern hemisphere
    slp  = slp.sel( shidx).load()
    topo = topo.sel(shidx).load()

    # adjust slp units
    slp = slp.where(False, slp / 100)
    slp.attrs['units'] = 'hPa'

    # adjust topography data
    topo = topo.where(False, topo / 9.8)            # geop. to geop. height
    topo = topo.sel({'time': topo['time'][[0]]})     # leave first timestep
    topo.attrs['units'] = 'm'

    # convert coordinates to lambert projection
    slp_xy  = convert_to_lambert(slp)
    topo_xy = convert_to_lambert(topo).squeeze()

# find cyclone centers
# ____________________

def slp_cyclonic_minima(slp_xy, datestr):
    '''
    info:
    parameters:
    returns:
    '''

    # initiate timestep subcontainer
    cyclones = []

    # process each latitude gridpoint
    for i in range(1, len(slp_xy['y'])-1):

        # process each longitude gridpoint
        for j in range(1, len(slp_xy['x'])-1):

            # 9 points box (point in index 4 is the center)
            slp_xy_box  = slp_xy[ i-1:i+2, j-1:j+2].data.flatten()
            topo_xy_box = topo_xy[i-1:i+2, j-1:j+2].data.flatten()

            # separate points
            slp_xy_boxcenter  = slp_xy_box[4]
            slp_xy_boxborders = np.delete(slp_xy_box, 4)
            topo_xy_boxcenter = topo_xy_box[4]

            # initial conditions to consider minimum as cyclone
            cond1 = ((slp_xy_boxcenter  <  slp_xy_boxborders).sum() == 8)
            cond2 = ( topo_xy_boxcenter <= topo_max)

            # check if center is minimum
            if cond1 and cond2:

                # add x-y coordinates to subcontainer
                cyclones.append([slp_xy[i, j]['y'].item(),
                                 slp_xy[i, j]['x'].item(),
                                 slp_xy_boxcenter])

    # process each point previously selected
    ipoint = 0
    while ipoint < len(cyclones):

        # separate point coordinates and values
        yp    = cyclones[ipoint][0]
        xp    = cyclones[ipoint][1]
        slp_p = cyclones[ipoint][2]

        # grids center in point
        XC = (slp_xy['X'] - xp).data
        YC = (slp_xy['Y'] - yp).data

        # mask that defines 1000 km area surrounding point
        mask = ((XC**2 + YC**2) < radius**2)

        # container for gradient values
        grad = []

        # container for border of search radius
        slp_circle_border = []

        # process circle border
        for j in range(len(slp_xy['x'])):

            # check if y-coordinate has, at least, one value to retrieve
            if (mask[j, :].sum() < 2):

                # skip y-coordinate
                continue

            # # border coordinates
            # yb = slp_xy_k['Y'][j, mask[j, :]][[0, -1]].data
            # xb = slp_xy_k['X'][j, mask[j, :]][[0, -1]].data
            #
            # # process both border points
            # for pb in zip(yb, xb):
            #
            #     # add border point to container
            #     slp_circle_border.append((yb, xp))

            # gradients in circle border
            grad.append((slp_xy[j, mask[j, :]][ 0] - slp_p).item())
            grad.append((slp_xy[j, mask[j, :]][-1] - slp_p).item())

        # mean gradient around center
        gradm = np.nanmean(grad)

        # next conditions to consider minimum as cyclone
        cond = (gradm >= grad_min)

        # check if center is minimum (if not, remove from container)
        if not cond:

            # remove point
            _ = cyclones.pop(ipoint)

        else:

            # add gradient to point container and border points of search radius
            cyclones[ipoint] += [gradm, datestr]

            # continue to next point
            ipoint += 1

    # output
    return cyclones


def slp_cyclonic_minima_wrapper(args):
    '''
    info:
    parameters:
    returns:
    '''

    # run main function
    return slp_cyclonic_minima(*args)


def xarray_time_iterator(da):
    '''
    info:
    parameters:
    returns:
    '''

    # create iterator
    for t in da['time']:

        # define date string
        datestr = t.dt.strftime('%Y-%m-%d').item()

        # generate iterator item
        yield (da.sel({'time': t}), datestr)


# load cyclones if possible
if (not calc) and (os.path.isfile(f'{dirout}/{fout}')):

    # open cyclones file
    handle   = open(f'{dirout}/{fout}', 'rb')
    cyclones = pickle.load(handle)
    handle.close()

# calculate cyclone center positions
else:

    # cyclone container
    cyclones = {}

    # slp minimum identifier
    sid = 0

    # create thread pool
    pool = mp.Pool(processes=25)

    # arguments for minimum selection process functioimport multiprocessing as mpn
    args = xarray_time_iterator(slp_xy)

    # compute processing
    results = pool.imap(slp_cyclonic_minima_wrapper, args)

    # create cyclones container
    cyclones = {}

    # fill container
    for k, cyclones_k in enumerate(results):

        # key to dictionary (date)
        datestr = slp_xy.isel({'time': k})['time'].dt.strftime('%Y-%m-%d').item()
        dateid  = slp_xy.isel({'time': k})['time'].dt.strftime(' %Y%m%d ').item()

        # add to dictionary
        cyclones[datestr] = cyclones_k

        # add cyclone identifier to points
        for point in cyclones[datestr]:

            point.append(sid)
            sid +=1

        # logging message
        indent = log_message(f'slp centers found in {datestr}')
        print(f'{indent}{len(cyclones_k)} (second selection) ')

    # save dictionary with slp minimum centers in file
    handle = open(f'{dirout}/{fout}', 'wb')
    pickle.dump(cyclones, handle, protocol=pickle.HIGHEST_PROTOCOL)
    handle.close()

# plot cyclone centers
# ____________________

# execute only if asked
if plot:

    # output fiename template
    output = f'piss_map_slpmin_{simid.lower()}_*.png'

    # remove previous images
    for f in glob(f'{dirimg}/{output}'): os.remove(f)

    # levels of slp (for maps)
    levels = np.arange(970, 1050+1, 5)

    # plot tracks over map
    for t in slp_xy['time']:

        # date identifiers
        datestr = t.dt.strftime('%Y-%m-%d').item()
        dateid  = t.dt.strftime( '%Y%m%d' ).item()

        # time indexer
        tidx = {'time': t}

        # create figure
        fig, ax = plt.subplots(1, 1, figsize=size, subplot_kw={'projection': proj})

        # plot field
        slp.sel(tidx).plot.contourf(ax=ax,
                                    transform=trans,
                                    levels=levels,
                                    extent='both',
                                    cmap='jet',
                                    cbar_kwargs={'location'   : 'right',
                                                 'orientation': 'vertical',
                                                 'drawedges'  : False,
                                                 'fraction'   : 0.1,
                                                 'shrink'     : 1,
                                                 'aspect'     : 30,
                                                 'pad'        : 0.00,
                                                 'anchor'     : (0.5, 1),
                                                 'panchor'    : (0.5, 0)})

        # plot minimum points
        for p in cyclones[datestr]:

            # transform to lat-lon coordinates
            plat, plon = convert_point_to_latlon(p)

            # draw point in map
            ax.plot(plon, plat, lw=1, transform=trans,
                    marker='o', mfc='yellow', mec='black', zorder=4)#, ms=5)

            # write code of slp minimum
            ax.text(plon, plat, p[-1], size=12, transform=trans, zorder=5, color='gold')

            # # extract border points of search radius
            # latb = [pb[0] for pb in p[5]]
            # lonb = [pb[1] for pb in p[5]]
            #
            # # draw border points of search radius
            # ax.plot(lonb, latb, lw=0, transform=trans,
            #         marker='.', ms=1, mec='snow', mfc='snow')

        # add coastlines
        ax.coastlines(color='black', linewidth=0.5)

        # add and adjust gridlines
        gridlines = ax.gridlines(linewidth=1,
                                 color='grey',
                                 alpha=0.25,
                                 ls='--')

        gridlines.top_labels    = True
        gridlines.bottom_labels = True
        gridlines.left_labels   = True
        gridlines.right_labels  = True

        # set labels
        ax.set_title(f'Southern Hemishphere, {datestr}')
        ax.set_xlabel('')
        ax.set_ylabel('')

        # remove box
        # ax.axis('off')

        # set coordinate limits
        ax.set_extent([0, 360, -90, -9], crs=trans)

        # save / show plot
        fig.savefig(f"{dirimg}/{output.replace('*', dateid)}", dpi=300, bbox_inches='tight')
        plt.close()

