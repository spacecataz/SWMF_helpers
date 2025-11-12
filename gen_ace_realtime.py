#!/usr/bin/env python3

'''
Create an auto-updating SWMF solar wind input file that automatically updates
based on real-time ACE data.

Data is fetched from https://services.swpc.noaa.gov
Bad data values are removed via linear interpolation or, if occurring at end
of data stream, via persistence.

All data is ballistically propagated foward in time. Values newer than the
previous last good data point are appended to the output file.
'''

import urllib
import datetime
import json
from time import sleep
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import medfilt
from matplotlib.dates import date2num

from spacepy.pybats import ImfInput
from spacepy.datamodel import dmarray


# Start by configuring the argparser:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-o", "--outfile", type=str, default='IMF.dat',
                    help='Name of output file to dump data. ' +
                    'Defaults to IMF.dat')
parser.add_argument("-i", "--initfile", type=str, default=None,
                    help='Set initial file. Instead of starting from ' +
                    'scratch, append new values to the intial file.')
parser.add_argument("-s", "--smoothing", default=5, type=int,
                    help="Velocity may be smoothed via median filtering. " +
                    "Set this argument to an integer window size to apply " +
                    "smoothing. Default is 5 points (5 minute smoothing).")
parser.add_argument("-v", "--verbose", default=False, action='store_true',
                    help="Turn on verbose output mode.")
parser.add_argument("-w", "--wait", default=300, type=int,
                    help="Set wait time between refreshing data in seconds.")
# delay between pulls.
args = parser.parse_args()


address_swe = 'https://services.swpc.noaa.gov/text/ace-swepam.txt'
address_mag = 'https://services.swpc.noaa.gov/text/ace-magnetometer.txt'

allvars = ['n', 'ux', 't', 'bx', 'by', 'bz']
units = {v: u for v, u in zip(allvars, ['cm-3', 'km/s', 'K'] + 3*['nT'])}

# Declare important constants
RE = 6371              # Earth radius in kmeters.
l1_dist = 1495980      # L1 distance in km.
kboltz = 1.380649E-23  # Boltzmann constant, J/K
mp = 1.67262192E-27    # Proton mass in Kg
bound_dist = 32 * RE   # Distance to BATS-R-US upstream boundary from Earth.
travel_dist = l1_dist - bound_dist  # Actual distance to propagate.


def parse_json(stream):
    '''
    Given a json data stream from real time solar wind data,
    convert to a dictionary of Numpy arrays.

    Argument `stream` can accept a `.read()` compatable stream, e.g., file
    or URLresponse.
    '''

    raw = json.load(stream)
    header = raw.pop(0)

    # Create dict of varname:colnumber
    if 'bx_gsm' in header:
        varnames = {'bx': 1, 'by': 2, 'bz': 3}
    elif 'density' in header:
        varnames = {'n': 1, 'ux': 2, 't': 3}

    # Create data container.
    nlines = len(raw)
    data = {}
    data['time'] = np.zeros(nlines, dtype=object)
    for v in varnames.keys():
        data[v] = np.zeros(nlines)

    # Parse!
    for i, r in enumerate(raw):
        # Get time:
        data['time'][i] = datetime.datetime.strptime(r[0],
                                                     '%Y-%m-%d %H:%M:%S.000')

        # Get other variables:
        for v in varnames:
            data[v][i] = r[varnames[v]]

    return data


def parse_ascii(stream):
    '''
    Given a data stream from Real Time ACE, convert to a dictionary of
    Numpy arrays.

    Argument `stream` can accept a `.read()` compatable stream, e.g., file
    or URLresponse.
    '''

    # Read first line, set file type based on that.
    line = stream.readline()
    if type(line) is bytes:
        line = line.decode("utf-8")
    # Create dict of varname:colnumber
    if 'mag' in line:
        varnames = {'bx': 7, 'by': 8, 'bz': 9}
    elif 'swe' in line:
        varnames = {'n': 7, 'ux': 8, 't': 9}
    else:
        raise ValueError('Not a valid data stream.')

    # Skip past header.
    while '#-------' not in line:
        line = stream.readline()
        if type(line) is bytes:
            line = line.decode("utf-8")

    # Grab the rest of the lines.
    lines = stream.readlines()
    nlines = len(lines)

    # Create data container.
    data = {}
    data['time'] = np.zeros(nlines, dtype=object)
    for v in varnames.keys():
        data[v] = np.zeros(nlines)

    # Parse!
    for i, l in enumerate(lines):
        if type(l) is bytes:
            parts = l.decode('utf-8').split()
        else:
            parts = l.split()

        # Get time:
        tstring = ' '.join(parts[:4])
        data['time'][i] = datetime.datetime.strptime(tstring, '%Y %m %d %H%M')

        # Get other variables:
        for v in varnames:
            data[v][i] = parts[varnames[v]]

    return data


def fetch_rtsw(raw_swe=None, raw_mag=None):
    '''
    Fetch the current RT solar wind results and convert to SWMF format.
    '''

    if raw_swe is None:
        raw_swe = parse_ascii(urllib.request.urlopen(address_swe))
    if raw_mag is None:
        raw_mag = parse_ascii(urllib.request.urlopen(address_mag))

    # Check times: are we consistent?
    if (raw_swe['time'].size != raw_mag['time'].size) or \
       (raw_swe['time'][0] != raw_mag['time'][0]) or \
       (raw_swe['time'][-1] != raw_mag['time'][-1]):
        raise ValueError('SWEPAM and MAG times do not match.')

    # Combine data sets.
    data = raw_swe
    for v in allvars[3:]:
        data[v] = raw_mag[v]

    # Remove bad data values.
    time = date2num(data['time'])
    for v in allvars:
        # Get all valid points.
        locgood = data[v] > -999

        # If no valid points, raise exception.
        # If no bad points, no need to fill.
        if locgood.sum() == 0:
            raise ValueError(f'No valid data for variable {v}.')
        elif locgood.sum() == data['time'].size:
            continue

        # Linearly interpolate over bad values, using "last good value"
        # for the final points.
        tfilt, vfilt = time[locgood], data[v][locgood]
        func = interp1d(tfilt, vfilt, fill_value=(vfilt[0], vfilt[-1]),
                        bounds_error=False)
        data[v] = func(time)

    # Convert velocity to right coords:
    data['ux'] *= -1

    # Create seconds-from-start time array:
    tsec = np.array([(t - data['time'][0]).total_seconds()
                     for t in data['time']])

    # Time shift:
    velsmooth = medfilt(data['ux'], args.smoothing)
    shift = travel_dist / velsmooth
    tshift = np.array([t1 - datetime.timedelta(seconds=t2) for t1, t2 in
                      zip(data['time'], shift)])

    # Ensure that any points that are "overtaken" (i.e., slow wind overcome by
    # fast wind) are removed. First, locate those points:
    keep = [0]
    discard = []
    lasttime = tshift[0]
    for i in range(1, data['time'].size):
        if tshift[i] > lasttime:
            keep.append(i)
            lasttime = tshift[i]
        else:
            discard.append(i)
    if args.verbose:
        print(f'Removing "overtaken" points {len(discard)} of {tsec.size}.')

    # Create new IMF object and populate with propagated values.
    # Use the information above to throw out overtaken points.
    tnow = datetime.datetime.now()
    imfout = ImfInput(args.outfile+f'_{tnow:%Y%m%d_%H%M%S}',
                      load=False, npoints=len(keep))
    for v in allvars:
        imfout[v] = dmarray(data[v][keep], {'units': units[v]})
    imfout['time'] = tshift[keep]
    imfout.attrs['header'].append('Source data: ACE Real Time\n')

    return imfout


def test_fetch():
    '''
    This function is for verifying the data fetch function.
    '''

    fname1 = 'data/ace-swepam.txt'
    fname2 = 'data/ace-magnetometer.txt'
    with open(fname1, 'r') as f1, open(fname2, 'r') as f2:
        raw_swe = parse_ascii(f1)
        raw_mag = parse_ascii(f2)
        imf = fetch_rtsw(raw_swe, raw_mag)

    return imf


# ## START MAIN SCRIPT
# Initialize a file.
if args.initfile:
    print('Opening initial file...')
    imf = ImfInput(args.initfile)
    imf.attrs['file'] = args.outfile
else:
    print('Fetching initial data...')
    imf = fetch_rtsw()
    print('Success! Waiting for updates...')
    sleep(args.wait)

print('Beginning main loop.')

# Watch for data!
while True:
    print(20*'-')
    # Wait a bi

    tnow = datetime.datetime.now()
    print(f'UPDATING DATA AT T={tnow}')

    # Get updated data:
    print('\tFetching data....')
    imf_new = fetch_rtsw()
    print('\tSuccess. Saving interim file.')
    imf_new.write()

    # Update main data file:
    loc = imf_new['time'] > imf['time'][-1]
    if loc.sum() == 0:
        print('\nNo new data values.')
        sleep(args.wait)
        continue
    print(f'\tAppending {loc.sum()} new values...')
    imf['time'] = np.append(imf['time'], imf_new['time'][loc])
    for v in imf.attrs['var']:
        imf[v] = np.append(imf[v], imf_new[v][loc])

    print('\tWriting updated file.')
    imf.write()
    sleep(args.wait)
