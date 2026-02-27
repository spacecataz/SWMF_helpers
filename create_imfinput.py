#!/usr/bin/env python3

'''
Create an solar wind/IMF input file for use in SWMF-Geospace.

Input data can be fetched automatically or created from existing files.
Data will be cleaned and propagated from L1 to +32RE (the typical Geospace
upstream boundary).

Possible sources include ACE, DSCOVR, and WIND data through:

- CDAWeb's HAPI server
- NASA CDF files
- NOAA's real-time data streams
- Existing SWMF input files (useful for propagating unpropagated input data).

'''

from os import path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datetime import datetime, timedelta

from spacepy.plot import style

import sys

installdir = path.dirname(__file__)
sys.path.append(installdir)

import sw_tools

style()

parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("source", type=str,
                    help='The source of the data. For auto-fetching, set to ' +
                    'ACE, DSCOVR, or WIND. For NASA CDFs or SWMF input ' +
                    'use the CDF file name (must end in .cdf or .dat, ' +
                    'respectively).')
parser.add_argument("-f2", "--file2", type=str, help="Specify the matching " +
                    "file for two-file sets (e.g., ACE SWE and MFI files). " +
                    "Note that the script will auto-search for the 2nd file " +
                    "if this argument is not set.")
parser.add_argument("-v", "--verbose", default=False, action='store_true',
                    help="Turn on verbose output mode.")
parser.add_argument('-t1', '--start', help='Start time to begin fetching in ' +
                    'YYYYMMDD format.', type=str, default=None)
parser.add_argument('-t2', '--end', help='End time for data fetching in ' +
                    'YYYYMMDD format. If not given, will default to ' +
                    'start + 1 day', type=str, default=None)
parser.add_argument("--debug", default=False, action='store_true',
                    help="Turn on debugging mode.")
parser.add_argument("-o", "--outfile", default='IMF', help="Set " +
                    "output file name.  Defaults to 'IMF.dat'")
parser.add_argument("-s", "--smoothing", default=1, type=int,
                    help="Velocity may be smoothed via median filtering. " +
                    "Set this argument to an integer window size to apply " +
                    "smoothing. Default is 1 point, or no smoothing.")
parser.add_argument("--tshift", "-t", type=float, default=-1.0,
                    help="Use a constant time shift, given in minutes, " +
                    "instead of ballistic propagation. Default is -1.0, " +
                    "which uses default ballistic propagation. A value of " +
                    "0 sets no propagation.")
args = parser.parse_args()

# Add postfix to file if not given:
args.outfile = args.outfile + '.dat'*(args.outfile[-4:] != '.dat')

# STEP 1: Fetch or open data:
if args.source.lower() in ('ace', 'dscover', 'wind'):
    # Handle start/stop time:
    if args.start is None:
        raise ValueError('Must supply start time when fetching data from web.')
    tstart = datetime.strptime(args.start, '%Y%m%d')
    if args.end is None:
        tend = tstart + timedelta(days=1)
    else:
        tend = datetime.strptime(args.end, '%Y%m%d')

    # Fetch data:
    raw = sw_tools.fetch_hapi(args.source.lower(), tstart, tend,
                              outname=args.outfile + '_l1raw')

elif args[-4] == '.cdf':
    pass
elif args[-4] == '.dat':
    pass
else:
    raise ValueError('Unrecognized data source.')

# STEP 2: Propagate:
sw_tools.l1_propagate(raw, outfile=args.outfile, shift=args.tshift,
                      smoothwin=args.smoothing, verbose=args.verbose)
