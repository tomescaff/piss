#!/home/cr2/agomez/.conda/envs/conda/bin/python
# -*- coding: utf-8 -*-
#
# ~/projects/piss/src/piss_cyclone_tracking.py
#
# info  : track cyclones in the souther hemisphere in piss simulations (cesm2), based on
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


def distances_between_points(p1, p2):
    '''
    info:
    parameters:
    returns:
    '''

    # define ellipsoid
    g = pyproj.Geod(ellps='WGS84')

    # convert point to latlon coordinates
    lat1, lon1 = convert_point_to_latlon(p1)
    lat2, lon2 = convert_point_to_latlon(p2)

    # calculate distance
    az12, az21, dist = g.inv(lon1, lat1, lon2, lat2)

    # output (minimum distance in km)
    return dist / 1000


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
    dak = da.isel({'time': kt}).squeeze()

    # # process each contour
    # while ((ilat < len(lat)) and (ilon < len(lon))): # <-- option 1, i think it won't work
    #
    #     # for now, do nothing
    #     None


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
    D_fmt     = re.compile(r'^D\d')             # search distance

    # retrieve variables
    simid = [arg for arg in args if simid_fmt.match(str(arg))]   # sim. ID
    D     = [arg for arg in args if D_fmt.match(    str(arg))]   # search distance
    calc  = True if 'calc'  in args else False
    plot  = True if 'plot'  in args else False

    # check arguments
    simid = 'lgm_100' if not simid else simid[0]
    D     = 1000      if not D     else int(D[0][1:])

    # output
    return (simid, calc, D, plot)


###########################
##### LOCAL VARIABLES #####
###########################


# get variables from input arguments
simid, calc, D, plot = get_variables(sys.argv)


############################
##### LOCAL PARAMETERS #####
############################


# directories
homedir = os.path.expanduser('~')               # home directory
dirin   = f'/mnt/cirrus/results/friquelme'      # cesm simulations
dirout  = f'{homedir}/projects/piss/data'       # data output
dirimg  = f'{homedir}/projects/piss/img'        # output for figures

# filenames for stored results
fin  = f'cyclones_{simid}.pkl'
fout = f'tracks_{simid}.pkl'

# variable to process
key = 'PSL'

# indexer to choose southern hemisphere (from -30° to the south)
shidx = {'lat': slice(-90, -5)}

# parameters needed in method
#D  = 1300   # threshold distance to track cyclones [km]
D0 = D    # threshold distance to track cyclones [first timesteps]

# figure parameters
size = (7.5, 6)


###################
##### RUNNING #####
###################


# first check if data should be loaded
if ('slp' not in locals()):

    # load simulation datasets
    ds = load_simulation(dirin, simid, 'h1', ['PSL' ])

    # separate variables
    slp  = ds['PSL']

    # extract temporal range (only for slp)
    slp  = slp.sel( {'time': slice(date_ini, date_end)})

    # extract southern hemisphere
    slp  = slp.sel(shidx).load()

    # adjust slp units
    slp = slp.where(False, slp / 100)
    slp.attrs['units'] = 'hPa'

    # convert coordinates to lambert projection
    slp_xy  = convert_to_lambert(slp)

# check if cyclones center positions file exists
if not os.path.isfile(f'{dirout}/{fin}'):

    # run cyclone identification script
    os.system(f'{sys.executable} piss_cylone_identification.py calc')

# open cyclones file
handle   = open(f'{dirout}/{fin}', 'rb')
cyclones = pickle.load(handle)
handle.close()

# track cyclone centers
# _____________________
# * point = [y[0], x[1], slp[2], grad_slp[3], date[4], sid[5]]

# cyclone code and id
cn  = 0
cid = f'cyclone{cn:03d}'

# logging message
indent = log_message('assembling tracks')

# load tracks if possible
if (not calc) and (os.path.isfile(f'{dirout}/{fout}')):

    # open cyclones file
    handle = open(f'{dirout}/{fout}', 'rb')
    tracks = pickle.load(handle)
    handle.close()

# calculate cyclone tracks
else:

    # tracks container
    tracks = {}

    # process each date
    for it in range(len(slp_xy['time'])-2):

        # date identifiers for initial step of track
        datestr0 = slp_xy['time'][it+0].dt.strftime('%Y-%m-%d').item()
        datestr1 = slp_xy['time'][it+1].dt.strftime('%Y-%m-%d').item()
        datestr2 = slp_xy['time'][it+2].dt.strftime('%Y-%m-%d').item()

        # process each point of timestep
        ip0 = 0
        while ip0 < len(cyclones[datestr0]):

            # extract point
            p0 = cyclones[datestr0][ip0]

            # create containers for chosen points
            trackpoints = []
            idx         = []

            # minimum distance neccesary to choose points for track
            dist_min = D

            # if p0[-1] == 85:
            #
            #         print(f'p0 [{p0[-1]}]: ({p0[0]}, {p0[1]})')
            #         print(f'p1 [{p1[-1]}]: ({p1[0]}, {p1[1]})')
            #         print(f'p2 [{p2[-1]}]: ({p2[0]}, {p2[1]})')
            #         print(f'p2_guess: ({p2_guess[0]}, {p2_guess[1]})')
            #         print(f'dist: {dist_guess} ({east})')
            #         input()

            # process points of next timestep
            for ip1, p1 in enumerate(cyclones[datestr1]):

                # distance between points (different timesteps)
                d01  = distance_between_points(p0, p1)
                east = east_or_not(p0, p1)

                # check if distance is less than initial threshold
                if (d01 < D0) and (east):

                    # guess next point
                    p2_guess = guess_point(p0, p1)

                    # process points of third timestep
                    for ip2, p2 in enumerate(cyclones[datestr2]):

                        # distance between point and guess (same timestep)
                        dist_guess = distance_between_points(p2, p2_guess)

                        # check if point is east to previous point
                        east = east_or_not(p1, p2)


                        # check if distance is less than threshold
                        if (dist_guess < dist_min) and (east):

                            # update minimum distance
                            dist_min = dist_guess

                            # replace points in containers
                            trackpoints = [ p0,  p1,  p2]
                            idx         = [ip0, ip1, ip2]

            # check if track was initialized
            if not trackpoints:

                # if not initialized, continue to next point
                ip0 += 1
                continue

            # remove points from original container
            _ = cyclones[datestr0].pop(idx[0])
            _ = cyclones[datestr1].pop(idx[1])
            _ = cyclones[datestr2].pop(idx[2])

            # continue adding points to track
            for jt in range(it+2, len(slp_xy['time'])-1):

                # date identifier
                datestr_next = slp_xy['time'][jt+1].dt.strftime('%Y-%m-%d').item()

                # extraction of last points in track
                plast = trackpoints[-1]
                pprev = trackpoints[-2]

                # guess for next point in track
                pnext_guess = guess_point(pprev, plast)

                # minimum distance neccesary to choose points for track
                dist_min = D

                # new flag to check if continuation of track was successful
                success = False

                # process each point in next timestep
                for ipnext, pnext in enumerate(cyclones[datestr_next]):

                    # distance between point and guess (same timestep)
                    dist_guess = distance_between_points(pnext, pnext_guess)

                    # check if point is east to last trackpoint
                    east = east_or_not(plast, pnext)

                    # print(f'pprev  : {pprev[0]}, {pprev[1]}')
                    # print(f'plast  : {plast[0]}, {plast[1]}')
                    # print(f'pnext  : {pnext[0]}, {pnext[1]}')
                    # print(f'pnext* : {pnext_guess[0]}, {pnext_guess[1]}')
                    # print(f'dist   : {dist_guess:.2f} ({east})')
                    # print('')
                    # input()

                    # check if distance is less than threshold
                    if (dist_guess < dist_min) and (east):

                        # update minimum distance
                        dist_min = dist_guess

                        # replace point in containers
                        newpoint = pnext
                        idx      = ipnext

                        # update continuation flag
                        success = True

                # check if new point needs to be added to track
                if not success:

                    # add final track to tracks container
                    tracks[cid] = trackpoints

                    # logging message
                    print(f'{indent}{datestr0}: {cid} ({len(trackpoints):02d})')

                    # update cyclone code and id
                    cn += 1
                    cid = f'c{cn:05d}'

                    # if track ended, continue to next point
                    break

                # add new point to track
                trackpoints.append(newpoint)

                # remove point from original container
                _ = cyclones[datestr_next].pop(idx)

    # save dictionary with cyclone tracks in file
    handle = open(f'{dirout}/{fout}', 'wb')
    pickle.dump(tracks, handle, protocol=pickle.HIGHEST_PROTOCOL)
    handle.close()

# plot cyclone tracks
# ____________________

# check if plot tracks
if plot:

    # output filename template
    output = f'piss_map_tracks_{simid.lower()}_*.png'

    # remove previous images
    for f in glob(f'{dirimg}/{output}'): os.remove(f)

    # levels of slp (for maps)
    levels = np.arange(970, 1050+1, 5)

    # plot tracks over map
    for cid in tracks.keys():

        # extract date strings
        date_ini = tracks[cid][ 0][-2]
        date_end = tracks[cid][-1][-2]

        # calculate mean slp between dates
        slpm = slp.sel({'time': slice(date_ini, date_end)}).mean(dim='time', keep_attrs=True)

        # coordinates of track
        latlon = [convert_point_to_latlon(point) for point in tracks[cid]]
        lat    = np.array([ll[0] for ll in latlon])
        lon    = np.array([ll[1] for ll in latlon])

        # fix range of longitude
        if (lon.max() - lon.min()) > 180:

            lon[lon > 180] = lon[lon > 180] - 360

        # create figure
        fig, ax = plt.subplots(1, 1, figsize=size, subplot_kw={'projection': proj})

        # plot field
        slpm.plot.contourf(ax=ax,
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

        # draw points in map
        ax.plot(lon,  lat, transform=trans,
                lw=0, ls='-', color='black',
                marker='o', mfc='yellow', mec='black', ms=5)

        # put days (from starting point) over
        for k in range(len(lat)):

            # draw number of day
            ax.text(lon[k], lat[k]+0.45, f'{k+1}', transform=trans, size=10, color='black')

        # add coastlines
        ax.coastlines(color='black', linewidth=0.5)

        # add and adjust gridlines
        gridlines = ax.gridlines(linewidth=1,
                                 color='grey',
                                 alpha=0.25,
                                 ls='--')

        gridlines.bottom_labels = True
        gridlines.left_labels   = True

        # set labels
        ax.set_title(f'Cyclone track {cid[-3:]}, {date_ini} to {date_end}')
        ax.set_xlabel('')
        ax.set_ylabel('')

        # remove box
        # ax.axis('off')

        # set coordinate limits
        ax.set_extent([0, 360, -90, -9], crs=trans)

        # save / show plot
        fig.savefig(f"{dirimg}/{output.replace('*', cid)}", dpi=300, bbox_inches='tight')
        plt.close()


























