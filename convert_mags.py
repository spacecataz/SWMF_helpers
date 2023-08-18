#!/usr/bin/env python3

'''
This script converts magnetometer output from the SWMF into the format required
by CCMC and by the SWPC validation suite.  A new file is made for each
magnetometer.

If the Supermag package is installed, the real lat and lon of each station
will be included int he output file. If not, a dummy value will be used.
https://github.com/spacecataz/supermag
'''

import os
import datetime as dt
from glob import glob
from argparse import ArgumentParser

# Default values:
magfile = None

parser = ArgumentParser(description=__doc__)

parser.add_argument("file", type=str, default='./',
                    help='A path to a single magnetometer output file or ' +
                    'directory containing a magnetometer output file.')
parser.add_argument('-o', '--outdir', type=str, default='./', help='Path to '
                    + 'location to place output files.  Defaults to PWD.  ' +
                    'If outdir does not exist, create it.')
parser.add_argument('-d', '--debug', action='store_true',
                    help='Turn on debug info.')
parser.add_argument('-c', '--comp', action='store_true', help='Add full '
                    + 'contribution breakdown (e.g., dBnMhd, dBnFac, etc.)')

# Parse arguments and stash magfile into convenience variable:
args = parser.parse_args()
magfile = args.file

# Check outdir; create as necessary.
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)
if args.debug:
    print(f'Saving result to {args.outdir}')

# If a directory is given, search for the mag file:
if os.path.isdir(magfile):
    found_files = glob(magfile+'/magnetometer*.mag')
    if not found_files:
        raise ValueError('No magnetometer file found in {}.'.format(magfile))
    magfile = found_files[0]

if args.debug:
    print(f'Converting magnetometer file {magfile}...')

# Open magnetometer file for reading:
infile = open(magfile, 'r')

# Read header to get list of variables and stations.
stats = infile.readline().split(':')[-1].split()
varnames = infile.readline().split()

if args.debug:
    print('{} Stations found:'.format(len(stats)))
    for i, s in enumerate(stats):
        print('\t#{:04d}=={:}'.format(i+1, s))

# Slurp rest of file.
if args.debug:
    print('Reading mag file. This can take a moment...')
rawlines = infile.readlines()
infile.close()

try:
    import supermag
    info = supermag.read_statinfo()
except ModuleNotFoundError:
    info = {}

if info and args.debug:
    print('Loading station info successful. Using real lat/lons.')
else:
    print('Station info not found. Using dummy lat/lons.')

# Create new files, write headers.
for i, s in enumerate(stats):
    # Get lat/lon for file (used to be fixed at lon=355.310, lat=55.6300)
    if s in info:
        lat, lon = info[s]['aacgmlat'], info[s]['aacgmlon']
    else:
        lat, lon = 99.9999, 999.9999

    if args.debug:
        print('Working on station {}...'.format(s))
    f = open(args.outdir+'/{}.txt'.format(s), 'w')
    f.write('# SWMF run: SWMF_SWPC\n')
    f.write('#SWMF run finished on {}\n'.format(dt.datetime.now().isoformat()))
    f.write('# North, East and vertical components of magnetic field\n')
    f.write('# computed from magnetosphere & ionosphere currents\n')
    f.write('# Station: {}\n'.format(s.lower()))
    f.write(f'# Position (MAG): lon=       {lon:7.3f} lat=      {lat:7.4f}\n')
    f.write('Year Month Day Hour Min Sec GeomagLat GeomagLon B_NorthGeomag')
    f.write(' B_EastGeomag B_DownGeomag')
    if args.comp:
        f.write('dBnMhd dBeMhd dBdMhd dBnFac dBeFac dBdFac dBnHal dBeHal ' +
                'dBdHal dBnPed dBePed dBdPed')
    f.write('\n[year] [month] [day] [hour] [min]')
    f.write(' [s] [deg] [deg] [nT] [nT] [nT]\n')

    # Set the number of variables to write out:
    index_stop = 27 if args.comp else 15

    # Loop through lines related to this magnetometer:
    for line in rawlines[i::len(stats)]:
        # Parse line and turn into floating-point values.
        parts = line.split()
        parts = [float(p) for p in parts]

        # Write time:
        f.write('{1:4.0f}{2:5.0f}{3:5.0f}{4:5.0f}{5:5.0f}{6:5.0f}'.format(
                *parts))
        # Write lat-lon:
        f.write('{:13.3f}{:13.3f}'.format(0.0, 0.0))
        # Write perturbation:
        for p in parts[12:index_stop]:
            f.write(f'{p:13.3f}')
        f.write('\n')

    f.close()
