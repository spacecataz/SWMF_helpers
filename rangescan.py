#!/usr/bin/env python3
'''
Range Scan will scan a set of SWMF output files and report the range of
values within. This helps with, for example, determining what color bar
limits to use when plotting contours over a time series of files.

'''

import re
import argparse

import numpy as np

from spacepy.pybats import rim, bats

# Build, populate, and execute argument parser:
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("files", nargs='+', help="Files to scan; should be a " +
                    "class of SWMF output (e.g., RIM it*.dat files, MHD 2D " +
                    "cuts containing the same variable set, etc.). Can be " +
                    "explicit files or a unix wildcard.")
parser.add_argument("--debug", "-d", action="store_true",
                    help="Turn on debug information.")

args = parser.parse_args()

knowntypes = ['RimIono', 'Bats2d', 'IdlFile', 'MagGrid']

# Use the first file to guess at the proper parser:
filetype = 'IdlFile'  # Default.
filename = args.files[0].split('/')[-1]
if re.search('it.*\.idl(\.gz){0,1}', filename):
    filetype = 'RimIono'
    readfile = rim.Iono
elif re.search('mhd.*\.out(s)*', filename):
    filetype = 'Bats2d'
    readfile = bats.Bats2d
elif re.search('mag_grid', filename):
    filetype = 'MagGrid'
    readfile = bats.MagGridFile
else:
    print('\tUnrecognized file type; defaulting to IdlFile.')
    readfile = bats.IdlFile

if args.debug:
    print(f"File type detected is {filetype}")

# Get number of files, start tracking "bad" files.
nFile, nBad = len(args.files), 0

# Get list of variables in file. Remove common non-vars.
data = readfile(args.files[0])
varnames = list(data.keys())
if 'grid' in varnames:
    varnames.remove('grid')

# Create information containers:
nPts = np.zeros(nFile)
nPts[0] = data[varnames[0]].size
max, min, mean = {}, {}, {}
for v in varnames:
    max[v], min[v], mean[v] = np.zeros(nFile), np.zeros(nFile), np.zeros(nFile)
    max[v][0] = data[v].max()
    min[v][0] = data[v].min()
    mean[v][0] = data[v].mean()

# Loop through remaining files, continue scanning file:
for i, f in enumerate(args.files[1:]):
    if args.debug:
        print(f"Loading file {f}...")
    try:
        data = readfile(f)
    except Exception:
        if args.debug:
            raise IOError(f'{f} could not be opened.')
        print(f"Failed to load {f}")
        # If file fails to load, populate with nans:
        nBad += 1
        for v in varnames:
            min[v][i+1], max[v][i+1], mean[v][i+1] = 3 * [np.nan]
        continue

    for v in varnames:
        min[v][i+1] = data[v].min()
        max[v][i+1] = data[v].max()
        mean[v][i+1] = data[v].mean()
    nPts[i+1] = data[v].size

# Write results to screen:
print(f"RESULTS: {filetype}-type files ({nFile} found)")
print("Variable (Units)\tMin\tMax\tMean")
print(50*"-")
for v in varnames:
    meannow = np.sum(np.ma.masked_invalid(mean[v]))/(nFile-nBad)
    maxnow = np.ma.masked_invalid(max[v]).max()
    minnow = np.ma.masked_invalid(min[v]).min()
    units = 'No Units'
    if 'units' in data[v].attrs:
        units = data[v].attrs['units']
    print(f"{v} ({units}):\t{minnow:.3E}\t{maxnow:.3E}\t{meannow:.3E}")
