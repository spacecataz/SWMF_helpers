#!/usr/bin/env python3
'''
Propagate solar wind parameters from L1 to +32RE upstream for
use in the Space Weather Modeling Framework. Either a ballistic method or
constant time shift can be used. See the `tshift` argument for selecting
the propagation method.

Given a CDF file that contains the requisite values (see below), the time
array is delayed by the time each solar wind "parcel" (point in the file)
would take to travel from the L1 position to the upstream BATS-R-US boundary.
The travel time is taken using the parcel's current Earthward flow speed
(Vx), yielding a dynamic time adjustment.

If spacecraft position is found in the file, it is used to set the distance
over which the signal will be propagated dynamically. If not found, a static
distance will be used (distance to L1 point).

If plasma and field data have different sampling frequencies, the plasma
sampling times will be used. Values are mapped via linear interpolation.

CURRENT ASSUMPTIONS (to change):
- GSM coordinates for B, V vectors.
- The spacecraft is located along the GSM X-axis.

This script works with ACE and WIND data obtained from CDAWeb.

ACE spacecraft CDFs come in pairs: a SWEPAM file with plasma information,
spacecraft position, etc. and an MFI file with magnetic field data. If you
specify either file as the input argument, this script will attempt to find
the matching file in the same directory. On CDAWeb, select `AC_H0_MFI` and
`AC_H0_SWE` datasets and obtain values in GSM coordinates.
ACE spacecraft CDFs must contain the following variables:
| Variable Name(s)      | Description                                         |
|-----------------------|-----------------------------------------------------|
| SC_pos_GSM (optional) | Distance from Earth (km) in GSM/GSE coordinates     |
| alpha_ratio           | Alpha to proton ratio (not currently used)          |
| BGSM                  | Vector magnetic field (nT) data in GSM coordinates  |
| V_GSM                 | Vector solar wind velocity (km/s) in GSM coordinates|
| Np                    | Proton number density (1/ccm)                       |
| Tpr                   | Plasma temperature (K)                              |
| Epoch                 | CDF-formatted universal time.                       |

WIND spacecraft CDFs are single file downloads.
WIND spacecraft CDFs must contain the following variables:
| Variable Name(s) | Description                                            |
|------------------|--------------------------------------------------------|
| XGSM (optional)  | Distance from Earth (Re) in GSM/GSE coordinates.       |
| BX, BY, BZ       | Vector magnetic field (nT) data in GSM coordinates.    |
| VX, VY, VZ       | Vector solar wind velocity (km/s) in GSM coordinates.  |
| Np               | Proton number density (1/ccm)                          |
| TEMP             | Plasma temperature (K)                                 |
| Epoch            | CDF-formatted universal time.                          |
'''

from glob import glob
from os import path
import re
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datetime import datetime, timedelta

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import medfilt
from matplotlib.lines import Line2D

from spacepy.pycdf import CDF
from spacepy.plot import style
from spacepy.datamodel import dmarray
from spacepy.pybats import ImfInput

style()

# Start by configuring the argparser:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("file", type=str,
                    help='An L1 dataset in CDF format that contains the ' +
                    'required variables.')
parser.add_argument("-v", "--verbose", default=False, action='store_true',
                    help="Turn on verbose output mode.")
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

# Declare important constants
RE = 6371  # Earth radius in kmeters.
l1_dist = 1495980  # L1 distance in km.
bound_dist = 32 * RE  # Distance to BATS-R-US upstream boundary from Earth.
swmf_vars = ['bx', 'by', 'bz', 'ux', 'uy', 'uz', 'n', 't']
units = {v: u for v, u in zip(swmf_vars, 3*['nT']+3*['km/s']+['cm-3', 'K'])}


def pair(time1, data, time2, **kwargs):
    '''
    Use linear interpolation to pair two timeseries of data.  Data set 1
    (data) with time t1 will be interpolated to match time set 2 (t2).
    The returned values, d3, will be data set 1 at time 2.
    No extrapolation will be done; t2's boundaries should encompass those of
    t1.

    This function will correctly handle masked functions such that masked
    values will not be considered.

    **kwargs** will be handed to scipy.interpolate.interp1d
    A common option is to set fill_value='extrapolate' to prevent
    bounds errors.
    '''

    from numpy import bool_
    from numpy.ma import MaskedArray
    from scipy.interpolate import interp1d
    from matplotlib.dates import date2num

    # Dates to floats:
    t1 = date2num(time1)
    t2 = date2num(time2)

    # Remove masked values (if given):
    if type(data) is MaskedArray:
        if type(data.mask) is not bool_:
            d = data[~data.mask]
            t1 = t1[~data.mask]
        else:
            d = data
    else:
        d = data
    func = interp1d(t1, d, fill_value='extrapolate', **kwargs)
    return func(t2)


def read_wind(fname):
    # cdf_vars = ['BX', 'BY', 'BZ', 'VX', 'VY', 'VZ', 'Np', 'TEMP']
    # Open solar wind data and load variables into a dictionary:
    obs = CDF(fname)
    raw = {'bx': obs['BX'][...], 'by': obs['BY'][...], 'bz': obs['BZ'][...],
           'ux': obs['VX'][...], 'uy': obs['VY'][...], 'uz': obs['VZ'][...],
           'n': obs['Np'][...], 't': obs['TEMP'][...],
           'time': obs['Epoch'][...]}

    # If position: convert to km!
    return raw


def read_ace(fname):
    '''
    Given 1 of 2 ACE data CDFs (one for SWEPAM, one for MAG), find the matching
    CDF and load the data. Map to SWMF variables.

    Return a dictionary where each key is an SWMF var name mapped to the
    corresponding value in the ACE data.
    '''

    # Get critical parts of file name:
    result = re.search('ac\_h\ds\_mfi\_(\d+)\_(\d+)(\_cdaweb)?\.cdf',
                       args.file)
    stime, etime, cdatag = result.groups()
    if args.verbose:
        print(f"Found start/end times in file name of {stime}, {etime}")

    # Set path of files.
    dirname = path.dirname(fname)
    if dirname == '':
        dirname = './'

    # Find partner files.
    if 'swe' in fname:
        swe = CDF(fname)
        # Build matching name
        basefile = f'ac_h?s_mfi_{stime}_{etime}{cdatag}.cdf'
        fname2 = glob(path.join(dirname, basefile))[0]
        mag = CDF(fname2)
    elif 'mfi' in fname:
        mag = CDF(fname)
        # Build matching name
        basefile = f'ac_h?s_swe_{stime}_{etime}{cdatag}.cdf'
        fname2 = glob(path.join(dirname, basefile))[0]
        swe = CDF(fname2)
    else:
        raise ValueError('Expected "mfi" or "swe" in CDF file name.')

    # Extract plasma parameters from SWEPAM
    raw = {'time': swe['Epoch'][:], 'n': swe['Np'][:],
           't': swe['Tpr'][:], 'ux': swe['V_GSM'][:, 0],
           'uy': swe['V_GSM'][:, 1], 'uz': swe['V_GSM'][:, 2], }

    # Optional values:
    if 'SC_pos_GSM' in swe:
        raw['pos'] = swe['SC_pos_GSM'][:, 0]
    if 'alpha_ratio' in swe:
        raw['alpha'] = swe['alpha_ratio'][...]

    # Pair high-time resolution mag data to SWEPAM time:
    t_mag = mag['Epoch'][...]
    for i, b in enumerate(['bx', 'by', 'bz']):
        raw[b] = pair(t_mag, mag['BGSM'][:, i], raw['time'])

    return raw


# Open CDF and determine if this is an ACE or WIND data file.
if 'mfi' in args.file or 'swe' in args.file:
    raw = read_ace(args.file)

# Get S/C distance. If not in file, use approximation.
if 'pos' in raw:
    print('S/C location found! Using dynamic location.')
    raw['X'] = raw['pos'] - bound_dist
else:
    if args.verbose:
        print('S/C location NOT found, using static L1 distance.')
    raw['X'] = l1_dist - bound_dist

# Create seconds-from-start time array:
if args.tshift >= 0:
    if args.verbose:
        print(f'Using a STATIC timeshift of {args.tshift} minutes.')
        tsec = np.zeros(raw['time'].size) + args.tshift * 60.0
else:
    if args.verbose:
        print('Using BALLISTIC propagation.')
    tsec = np.array([(t - raw['time'][0]).total_seconds()
                     for t in raw['time']])

# Interpolate over bad data:
if args.verbose:
    print('Removing bad data:')
for v in swmf_vars + ['X']:
    # Find bad values:
    loc = ~np.isfinite(raw[v]) | (raw[v] <= -1E31)
    if args.verbose:
        print(f'\t{v} has {loc.sum()} bad values.')
    if loc.sum() == 0:
        continue
    interp = interp1d(tsec[~loc], raw[v][~loc], fill_value="extrapolate")
    raw[v][loc] = interp(tsec[loc])

# Apply velocity smoothing as required
if args.verbose:
    print(f'Applying smoothing using a {args.smoothing} window size.')
velsmooth = medfilt(raw['ux'], args.smoothing)

# Shift time: distance/velocity = timeshift (negative in GSM coords)
shift = raw['X']/velsmooth  # Time shift per point.
tshift = np.array([t1 - timedelta(seconds=t2) for t1, t2 in
                   zip(raw['time'], shift)])

# Ensure that any points that are "overtaken" (i.e., slow wind overcome by
# fast wind) are removed. First, locate those points:
keep = [0]
discard = []
lasttime = tshift[0]
for i in range(1, raw['time'].size):
    if tshift[i] > lasttime:
        keep.append(i)
        lasttime = tshift[i]
    else:
        discard.append(i)
if args.verbose:
    print(f'Removing "overtaken" points {len(discard)} of {tsec.size} total.')

# Create new IMF object and populate with propagated values.
# Use the information above to throw out overtaken points.
imfout = ImfInput(args.outfile+'.dat', load=False, npoints=len(keep))
for v in swmf_vars:
    imfout[v] = dmarray(raw[v][keep], {'units': units[v]})
imfout['time'] = tshift[keep]
imfout.attrs['header'].append(f'Source data: {args.file}')
imfout.attrs['header'].append('Ballistically propagted from L1 to upstream ' +
                              'BATS-R-US boundary')
imfout.attrs['header'].append(f'File created on {datetime.now()}')
imfout.write()

# Plot!
fig = imfout.quicklook(['by', 'bz', 'n', 't', 'ux'])
plotvars = ['by', 'bz', 'n', 't', 'ux']
for ax, v in zip(fig.axes, plotvars):
    c = ax.get_lines()[0].get_color()
    ax.plot(raw['time'], raw[v], '--', c=c, alpha=.5)
    ax.plot(raw['time'][...][discard], raw[v][...][discard],
            '.', c='crimson', alpha=.5)
l1 = Line2D([], [], color='gray', lw=4,
            label='Timeshifted Values')
l2 = Line2D([], [], color='gray', alpha=.5, linestyle='--', lw=4,
            label='Original Values')
l3 = Line2D([], [], marker='.', mfc='crimson', linewidth=0, mec='crimson',
            markersize=10, label='Removed Points')
fig.legend(handles=[l1, l2, l3], loc='upper center', ncol=3)
fig.subplots_adjust(top=.933)
fig.savefig(args.outfile + '_prop_info.png')
