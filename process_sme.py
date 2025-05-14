#!/usr/bin/env python3

'''
Given a SWMF magnetometer grid file, calculate SuperMag-like indexes and save
to file. If the --fetch flag is set, fetch matching indexes from SuperMag
and save along side SWMF values.
'''

import datetime as dt
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as np
from spacepy.datamodel import dmarray, SpaceData
from spacepy.pybats import bats
# import supermag_api as smapi


parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('magfile', type=str, help='Name of the magnetometer ' +
                    'grid file to process in to SM indexes.')

# Handle arguments:
args = parser.parse_args()

# Open our mag grid file, get key information:
mags = bats.MagGridFile(args.magfile)
nframe = mags.attrs['nframe']  # number of iterations for loop
shape = mags["dBn"].shape      # size of index used
# Get start/stop time.
tstart, tend = mags.attrs['time_range'][0], mags.attrs['time_range'][1]

# Create object to hold resulting data.
data = SpaceData()
data.attrs['descrip'] = 'Supermag-like indexes from SWMF MagGrid files.'

# Add time:
data['time'] = mags.attrs['times']

# Create empty arrays with units:
varnames = ['SWMFL', 'SWMFU', 'lat_L', 'lat_U', 'lon_L', 'lon_U']
for v in varnames:
    data[v] = dmarray(np.zeros(nframe))
    if 'SWMF' in v:
        data[v].attrs['units'] = 'nT'
    else:
        data[v].attrs['units'] = 'deg'

for i in range(nframe):
    mags.switch_frame(i)

    # Get location of max and min dBn:
    dbn = mags['dBn'].flatten()
    ind_min, ind_max = np.argmin(dbn), np.argmax(dbn)

    # Turn 1D index into 2D indices:
    ind_min = np.unravel_index(ind_min, shape)
    ind_max = np.unravel_index(ind_max, shape)

    # Save variables associated with max/min dBn:
    data['SWMFL'][i] = mags['dBn'][ind_min]
    data['SWMFU'][i] = mags['dBn'][ind_max]

    data['lon_L'][i] = mags['Lon'][ind_min[0]]
    data['lon_U'][i] = mags['Lon'][ind_max[0]]
    data['lat_L'][i] = mags['Lat'][ind_min[1]]
    data['lat_U'][i] = mags['Lat'][ind_max[1]]