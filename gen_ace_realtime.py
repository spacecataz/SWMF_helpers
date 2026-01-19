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

import datetime
import json
from time import sleep
from urllib.request import urlopen
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
parser.add_argument("-d", "--datasrc", default='rtsw', type=str,
                    help="Set data source: ace or rtsw.")
parser.add_argument("-v", "--verbose", default=False, action='store_true',
                    help="Turn on verbose output mode.")
parser.add_argument("-w", "--wait", default=300, type=int,
                    help="Set wait time between refreshing data in seconds.")
parser.add_argument("-r", "--refresh", default=20, type=int,
                    help="Set waiting time before refreshing real time data " +
                    "stream when an error is returned.")
# delay between pulls.
args = parser.parse_args()


address_swe = 'https://services.swpc.noaa.gov/text/ace-swepam.txt'
address_mag = 'https://services.swpc.noaa.gov/text/ace-magnetometer.txt'

address_rt = 'https://services.swpc.noaa.gov/products/solar-wind/'
# mag-6-hour.json

allvars = ['n', 'ux', 't', 'bx', 'by', 'bz']
units = {v: u for v, u in zip(allvars, ['cm-3', 'km/s', 'K'] + 3*['nT'])}

# Declare important constants
badflag = -999.0       # Bad data flag.
RE = 6371              # Earth radius in kmeters.
l1_dist = 1495980      # L1 distance in km.
kboltz = 1.380649E-23  # Boltzmann constant, J/K
mp = 1.67262192E-27    # Proton mass in Kg
bound_dist = 32 * RE   # Distance to BATS-R-US upstream boundary from Earth.
travel_dist = l1_dist - bound_dist  # Actual distance to propagate.


def unify_time(time1, time2):
    '''
    Given two timeseries, combine all unique points into a single array.
    '''

    time1 = time1.tolist()
    time2 = time2.tolist()

    for t in time2:
        if t not in time1:
            time1.append(t)

    time1.sort()
    return np.array(time1)


def pair(time1, data, time2, **kwargs):
    '''
    Use linear interpolation to pair two timeseries of data.  Data set 1
    (data) with time t1 will be interpolated to match time set 2 (t2).
    The returned values, d3, will be data set 1 at time 2.
    No extrapolation will be done; t2's boundaries should encompass those of
    t1.

    Bad data values will be removed prior to interpolation.

    **kwargs** will be handed to scipy.interpolate.interp1d
    A common option is to set fill_value='extrapolate' to prevent
    bounds errors.
    '''

    # Dates to floats:
    t1 = date2num(time1)
    t2 = date2num(time2)

    # Search for bad data values and remove:
    loc = ~np.isfinite(data) | (data <= badflag)

    # Trim down.
    t1, data = t1[~loc], data[~loc]
    if data.size == 0:
        raise ValueError('No valid data points in this time range.')

    # Create interpolator function
    func = interp1d(t1, data, fill_value=(data[0], data[-1]),
                    bounds_error=False, **kwargs)

    # Interpolate and return.
    return func(t2)


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


def fetch_rtsw(raw_swe=None, raw_mag=None, source='rtsw', duration=2):
    '''
    Fetch the current RT solar wind results and convert to SWMF format.
    Can fetch either ACE real time or combined ACE/DSCOVR

    Parameters
    ----------
    raw_swe, raw_mag : Data stream-like objects, defaults to None
        Set the data source. This is useful for testing against static file
        inputs.
    source : str, defaults to 'rtsw'
        Set data source, either 'ace' or 'rtsw' for ACE or combined
        ACE/DSCOVR Real Time Solar Wind data.
    duration : int, defaults to 2
        Set the amount of rtsw data to download, in hours. Can be 2 or 6.
    '''

    # Set address and parse strategy
    if source == 'ace':
        url_mag, url_pls = address_mag, address_swe
        parse = parse_ascii
    elif source == 'rtsw':
        url_mag = address_rt + f'mag-{duration}-hour.json'
        url_pls = address_rt + f'plasma-{duration}-hour.json'
        parse = parse_json

    # Data stream not provided, create:
    if raw_swe is None:
        raw_swe = parse(urlopen(url_pls))
    if raw_mag is None:
        raw_mag = parse(urlopen(url_mag))

    # We want no extrapolation into the future of any variable.
    # Remove trailing bad data values by finding the last good point for all
    # existing variables.
    # First, plasma data:
    locgood = True + np.zeros(raw_swe['time'].size, dtype=bool)
    for v in allvars[:3]:
        locgood = (locgood) & (raw_swe[v] > badflag)
    lastgood_p = raw_swe['time'][locgood][-1]
    # Then, mag data:
    locgood = True + np.zeros(raw_mag['time'].size, dtype=bool)
    for v in allvars[3:]:
        locgood = (locgood) & (raw_mag[v] > badflag)
    lastgood_m = raw_mag['time'][locgood][-1]

    # Find common last-good point:
    lastgood = min(lastgood_m, lastgood_p)
    print(f'\tLast valid data point is {lastgood}')

    # Trim down to last good data point
    loc = raw_swe['time'] <= lastgood
    for v in ['time'] + allvars[:3]:
        raw_swe[v] = raw_swe[v][loc]
    loc = raw_mag['time'] <= lastgood
    for v in ['time'] + allvars[3:]:
        raw_mag[v] = raw_mag[v][loc]

    # Get unified time:
    t_swe, t_mag = raw_swe['time'], raw_mag['time']
    time = unify_time(t_swe, t_mag)

    # Finally, pair the two data sets into one time series with no gaps.
    data = {'time': time}
    for v in allvars[:3]:
        data[v] = pair(raw_swe['time'], raw_swe[v], time)
    for v in allvars[3:]:
        data[v] = pair(raw_mag['time'], raw_mag[v], time)

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
        imf1 = fetch_rtsw(raw_swe, raw_mag)

    fname1 = 'data/plasma-6-hour.json'
    fname2 = 'data/mag-6-hour.json'
    with open(fname1, 'r') as f1, open(fname2, 'r') as f2:
        raw_swe = parse_json(f1)
        raw_mag = parse_json(f2)
        imf2 = fetch_rtsw(raw_swe, raw_mag)

    return imf1, imf2


# ## START MAIN SCRIPT
# Initialize a file.
if args.initfile:
    print(f'Opening initial file: {args.initfile}')
    imf = ImfInput(args.initfile)
    imf.attrs['file'] = args.outfile
else:
    print('Fetching initial data...')
    imf = fetch_rtsw(source=args.datasrc, duration=6)
    imf.attrs['file'] = args.outfile
    print(f'Writing initial values to {imf.attrs["file"]}')
    imf.write()
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
    try:
        imf_new = fetch_rtsw(source=args.datasrc, duration=2)
    except:
        print('\tWeb exception detected. Retrying...')
        sleep(args.refresh)
        continue
    print(f'\tSuccess. Saving interim file to {imf_new.attrs["file"]}')
    imf_new.write()

    # Update main data file:
    loc = imf_new['time'] > imf['time'][-1]
    if loc.sum() == 0:
        print('\tNo new data values.')
        sleep(args.wait)
        continue
    print(f'\tAppending {loc.sum()} new values to {imf.attrs["file"]}')
    imf['time'] = np.append(imf['time'], imf_new['time'][loc])
    for v in imf.attrs['var']:
        imf[v] = np.append(imf[v], imf_new[v][loc])

    print(f'\tWriting updated file ({imf.attrs["file"]}).')
    imf.write()
    sleep(args.wait)
