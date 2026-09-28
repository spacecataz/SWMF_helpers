#!/usr/bin/env python3

'''
Given a SWMF magnetometer grid .outs file, walk through each frame and
calculate the maximum field perturbation and dB/dt as a function of latitude
and longitude.

Currently, two values are saved: delta-B_H (the total horizontal component
of the surface magnetic field) and dB_H/dt (the magnitude of the rate of
change of the horizontal field). Note the latter uses the definition used in
the 2013 CCMC-SWPC validation challenge:

dB_h/dt = sqrt((dB_N/dt)^2 + (dB_E/dt)^2)

Additionally, the run time (in seconds from start) of the maximum value
is saved on the grid.

Results are saved as a dictionary stashed into a  Python pickle and can be
loaded via:

```
from pickle import load

with open('someresults.pkl', 'rb') as f:
    data = load(f)
```
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pickle import dump

import numpy as np
from spacepy.pybats.bats import MagGridFile

parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('fname', type=str, help='File path/name of mag grid ' +
                    'file to examine.')
parser.add_argument('-o', '--outfile', type=str, default='maggrid_max.pkl',
                    help="Name of output data file (should end in .pkl).")

# Handle arguments:
args = parser.parse_args()

# Open the file...
print(f"Opening mag grid file: {args.fname}")
mag = MagGridFile(args.fname)
mag.switch_frame(0)
mag.calc_h()

# Get some basic information, create data object:
nframe = mag.attrs['nframe']
nlon, nlat = mag['dBh'].shape
data = {'dbh': np.zeros([nlon, nlat]), 'dbt': np.zeros([nlon, nlat])}

# Save time and run time:
data['time_h'] = np.zeros([nlon, nlat])
data['time_t'] = np.zeros([nlon, nlat])

# Save coordinates:
data['lon'] = mag['Lon']
data['lat'] = mag['Lat']

# Arrays to hold data for comparison:
dbt_now = np.zeros([nlon, nlat])
dbn_last, dbe_last = np.zeros([nlon, nlat]), np.zeros([nlon, nlat])
dbn_last, dbe_last = mag['dBn'], mag['dBe']

# ...and away we go!
for i in range(1, nframe):
    print(f"Calculating maximums... {i/nframe:6.2%}", end='\r', flush=True)
    # Switch time to current:
    mag.switch_frame(i)

    # Calc dbdt:
    dt = mag.attrs['runtimes'][i] - mag.attrs['runtimes'][i-1]
    dn, de = mag['dBn'] - dbn_last, mag['dBe'] - dbe_last
    dbt_now = np.sqrt((dn/dt)**2 + (de/dt)**2)

    # Stash times for biggest values:
    data['time_h'][mag['dBh'] > data['dbh']] = mag.attrs['runtime']
    data['time_t'][dbt_now > data['dbt']] = mag.attrs['runtime']

    # Stash biggest values:
    data['dbh'] = np.fmax(data['dbh'], mag['dBh'])
    data['dbt'] = np.fmax(data['dbt'], dbt_now)

    # Update "last" values:
    dbn_last, dbe_last = mag['dBn'], mag['dBe']

# Save our data values:
with open(args.outfile, 'wb') as f:
    dump(data, f)