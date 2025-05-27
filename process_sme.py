#!/usr/bin/env python3

'''
Given a SWMF magnetometer grid file, calculate SuperMag-like indexes and save
to file. If the --fetch flag is set, fetch matching indexes from SuperMag
and save along side SWMF values.

The result is written to JSON-headed ASCII; the default file name is
'sm_indexes.txt'. The file can be read using spacepy as follows:

```
from spacepy import datamodel
data = datamodel.readJSONheadedASCII('sm_indexes.txt')
```

NOTE: The time is in a number format, to use dates use matplotlib as follows:

```
from matplotlib.dates import num2date
data['time'] = num2date(data['time'])
```

'''

import re
import datetime as dt
from matplotlib.dates import date2num  # num2date
import warnings
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from scipy import interpolate
import numpy as np
from spacepy.datamodel import dmarray, SpaceData
from spacepy.pybats import bats
import supermag_api as smapi


parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('magfile', type=str, help='Name of the magnetometer ' +
                    'grid file to process in to SM indexes.')
parser.add_argument("-o", "--outfile", default='sm_indexes.txt', help="Set " +
                    "output file name.  Defaults to 'sm_indexes.txt'")

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

# Add attributes to data object based on magnetometer grid input file:
data.attrs['header'] = mags.attrs['header'].strip()
# stash max/min lat in file as an attribute!

# Extract coordinate system for magnetometer location:
match = re.search('\((\w{3})\)', mags.attrs['header'])
data.attrs['coord'] = match.groups()[0]

# Warn if we are not in SMG coordinates.
# NOTE: For future use, be prepared to transform GEO or MAG to SMG.
if data.attrs['coord'] != 'SMG':
    warnings.warn('Input file not in SMG coordinates ' +
                  f'(coords={data.attrs["coord"]})')

# CHECK MAX/MIN LAT HERE! WARN AS NECESSARY!
# ??? Are degrees used
if mags['Lat'].max() != 80. or mags['Lat'].min() != 40.:
    warnings.warn('Latitude is outside of standard ' +
                  f"expected range = [40.0{chr(176)} - 80.0{chr(176)}]" +
                  f' actual range = [{mags["Lat"].min()}{chr(176)} - ' +
                  f'{mags["Lat"].max()}{chr(176)}]')
data.attrs['lat range'] = (f'{mags["Lat"].min()}{chr(176)} - ' +
                           f'{mags["Lat"].max()}{chr(176)}')

# Add time:
data['time'] = mags.attrs['times']

# Create empty arrays with units:
varnames = ['SWMFL', 'SWMFU', 'mlat_L', 'mlat_U', 'lon_L', 'lon_U']
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
    data['mlat_L'][i] = mags['Lat'][ind_min[1]]
    data['mlat_U'][i] = mags['Lat'][ind_max[1]]

# Convert SM longitude to MLT:
data['mlt_L'] = dmarray(data['lon_L'] / 15. + 12., {'units': 'Hours'})
data['mlt_U'] = dmarray(data['lon_U'] / 15. + 12., {'units': 'Hours'})

# Adjust MLT to account for the earth being a sphere
for i in range(nframe):
    if data['mlt_L'][i] > 24.:
        data['mlt_L'][i] = data['mlt_L'][i] - 24.

    if data['mlt_U'][i] > 24.:
        data['mlt_U'][i] = data['mlt_U'][i] - 24.

# Grab data from supermag. Extend range to ensure continuous interpolation.
SuperMag = smapi.fetch_index(tstart - dt.timedelta(minutes=+2),
                             tend + dt.timedelta(minutes=+2), 'amland')

# for testing as to not call the SuperMAG site everytime
SuperMag.toJSONheadedASCII('SuperMAG_' + args.outfile)

# Interpolate Data and Store in `data`!
time_real = np.linspace(date2num(data['time']).min(),
                        date2num(data['time']).max(), len(data['time']))
# Variables we want to include from SuperMAG:
SMvars = ['SMU', 'SMUmlat', 'SMUmlt', 'SML', 'SMLmlat', 'SMLmlt']
# Interpolate each in turn.
for v in SMvars:
    SM2SWMF = interpolate.interp1d(date2num(SuperMag['time']), SuperMag[v],
                                   kind='linear')
    data[v] = SM2SWMF(time_real)

# Add attributes from data
data.attrs['dates'] = (f'{data["time"].min()}' + ' - ' +
                       f'{data["time"].max()}')

# Fix times
data['time'] = date2num(data['time'])

# Write output file to disk.
data.toJSONheadedASCII(args.outfile)
