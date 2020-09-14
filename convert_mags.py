#!/usr/bin/env python

'''
This script converts magnetometer output from the SWMF into the format required
by CCMC and by the SWPC validation suite.  A new file is made for each
magnetometer.

Issues: need to incorporate station lat/lon into output.
'''

# Ensure that we are using a version of python >= 2.7.
import sys
if sys.version_info < (2,7):
    print('ERROR: Python version must be >= 2.7')
    print('Current version: '+sys.version)
    exit()
        

import os
import datetime as dt
from glob     import glob
from argparse import ArgumentParser

# Default values:
magfile = None

parser = ArgumentParser(description=__doc__)

parser.add_argument('-m', '--magfile', type=str, default='./',
                    help='A path to a single magnetometer output file or '
                    +' directory containing an output file.  Defaults to PWD.')
parser.add_argument('-o', '--outdir', type=str, default='./', help='Path to '
                    +'location to place output files.  Defaults to PWD.  '+
                    'If outdir does not exist, create it.')
parser.add_argument('-d', '--debug', action='store_true',
                    help='Turn on debug info.')
# Parse arguments and stash magfile into convenience variable:
args   = parser.parse_args()
magfile= args.magfile

# Check outdir; create as necessary.
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

# If a directory is given, search for the mag file:
if os.path.isdir(magfile):
    found_files = glob(magfile+'/magnetometer*.mag')
    if not found_files:
        raise(ValueError('No magnetometer file found in {}.'.format(magfile)))
    magfile = found_files[0]

if args.debug:
    print('Working on file {}'.format(magfile))
    
# Open magnetometer file for reading:
infile = open(magfile, 'r')

# Read header to get list of variables and stations.
stats    = infile.readline().split(':')[-1].split()
varnames = infile.readline().split()

if args.debug:
    print('{} Stations found:'.format(len(stats)))
    for i, s in enumerate(stats):
        print('\t#{:04d}=={:}'.format(i+1,s))
          

# Slurp rest of file.
rawlines = infile.readlines()
infile.close()

# Create new files, write headers.
files = []
for i,s in enumerate(stats):
    if args.debug:
        print('Working on station {}...'.format(s))
    f = open(args.outdir+'/{}.txt'.format(s), 'w')
    f.write('# SWMF run: SWMF_SWPC\n')
    f.write('#SWMF run finished on {}\n'.format(dt.datetime.now().isoformat()))
    f.write('# North, East and vertical components of magnetic field\n')
    f.write('# computed from magnetosphere & ionosphere currents\n')
    f.write('# Station: {}\n'.format(s.lower()))
    f.write('# Position (MAG): lon=       355.310 lat=      55.6300\n')
    f.write('Year Month Day Hour Min Sec GeomagLat GeomagLon B_NorthGeomag')
    f.write(' B_EastGeomag B_DownGeomag\n[year] [month] [day] [hour] [min]')
    f.write(' [s] [deg] [deg] [nT] [nT] [nT]\n')

    # Loop through lines related to this magnetometer:
    for l in rawlines[i::len(stats)]:
        # Parse line and turn into floating-point values.
        parts = l.split()
        parts = [float(p) for p in parts]

        # Write time:
        f.write('{1:4.0f}{2:5.0f}{3:5.0f}{4:5.0f}{5:5.0f}{6:5.0f}'.format(*parts))
        # Write lat-lon:
        f.write('{:13.3f}{:13.3f}'.format(0.0, 0.0))
        # Write perturbation:
        f.write('{12:13.3f}{13:13.3f}{14:13.3f}\n'.format(*parts))

    f.close()
