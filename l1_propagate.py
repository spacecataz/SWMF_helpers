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

This script works with ACE, Wind, and DSCOVR data obtained from CDAWeb.

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

WIND spacecraft CDFs are multi-file downloads; look for **wi_h1s_swe** and
**wi_h0s_mfi** file types on CDAWeb.
| Variable Name(s) | Description                                            |
|------------------|--------------------------------------------------------|
| PGSM (optional)  | Distance from Earth (Re) in GSM/GSE coordinates.       |
| BGSM             | Vector magnetic field (nT) data in GSM coordinates.    |
| Proton_VX_moment | Vector solar wind velocity (km/s) in GSE coordinates.  |
| Proton_VY_moment | ... in the Y direction                                 |
| Proton_VZ_moment | ... in the Z direction                                 |
| Np               | Proton number density (1/ccm)                          |
| Proton_W_moment  | Plasma thermal velocity, converted to temperature (K)  |
| Epoch            | CDF-formatted universal time.                          |

It is possible to obtain single-file formats for WIND; such CDFs must contain
the following variables:
| Variable Name(s) | Description                                            |
|------------------|--------------------------------------------------------|
| XGSM (optional)  | Distance from Earth (Re) in GSM/GSE coordinates.       |
| BX, BY, BZ       | Vector magnetic field (nT) data in GSM coordinates.    |
| VX, VY, VZ       | Vector solar wind velocity (km/s) in GSM coordinates.  |
| Np               | Proton number density (1/ccm)                          |
| TEMP             | Plasma temperature (K)                                 |
| Epoch            | CDF-formatted universal time.                          |

DSCOVR data come in two CDFs: one for the Faraday cup instrument (plasma) and
one for the magnetic field instrument. CDAweb files should begin with
dscovr_h0s_mag or dscovr_h1s_fc.
| Variable Name(s) | Description                                            |
|------------------|--------------------------------------------------------|
| B1GSE            | Vector magnetic field (nT) data in GSE coordinates.    |
| V_GSE            | Vector solar wind velocity (km/s) in GSE coordinates.  |
| Np               | Proton number density (1/ccm)                          |
| THERMAL_TEMP     | Plasma thermal temperature                             |
| Epoch/Epoch1     | CDF-formatted universal time in the fc/mag files.      |
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
parser.add_argument("-f2", "--file2", type=str, help="Specify the matching " +
                    "file for two-file sets (e.g., ACE SWE and MFI files). " +
                    "Note that the script will auto-search for the 2nd file " +
                    "if this argument is not set.")
parser.add_argument("-v", "--verbose", default=False, action='store_true',
                    help="Turn on verbose output mode.")
# parser.add_argument("-r", "--rotate", default=True,
#                     help="Rotate coordinates from GSE to GSM (default) as " +
#                     "applicable. If set to false, GSE coordinates are " +
#                     "used and noted in the file header.")
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
RE = 6371              # Earth radius in kmeters.
l1_dist = 1495980      # L1 distance in km.
kboltz = 1.380649E-23  # Boltzmann constant, J/K
mp = 1.67262192E-27    # Proton mass in Kg
bound_dist = 32 * RE   # Distance to BATS-R-US upstream boundary from Earth.

# Set var names and units.
swmf_vars = ['bx', 'by', 'bz', 'ux', 'uy', 'uz', 'n', 't']
units = {v: u for v, u in zip(swmf_vars, 3*['nT']+3*['km/s']+['cm-3', 'K'])}


def gse_to_gsm(x, y, z, time):
    '''
    Use spacepy's Coord module to rotate from GSE to GSM.

    Parameters
    ----------
    x, y, z : Numpy array-like
        The GSE X, Y, and Z values of the timeseries to rotate.
    time : array of Datetimes
        The time corresponding to the xyz timeseries.

    Returns
    -------
    x, y, z : Numpy arrays
        The rotated values now in GSM coordinates.
    '''

    from spacepy.coordinates import Coords
    from spacepy.time import Ticktock

    # Convert time to ticktocks:
    ticks = Ticktock(time, 'ISO')

    # Rearrange values to correct shape:
    xyz = np.array([x, y, z]).transpose()

    # Rotate:
    gse_vals = Coords(xyz, 'GSE', 'car', ticks=ticks)
    gsm_vals = gse_vals.convert('GSM', 'car')

    return gsm_vals.x, gsm_vals.y, gsm_vals.z


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


def pair(time1, data, time2, varname=None, **kwargs):
    '''
    Use linear interpolation to pair two timeseries of data.  Data set 1
    (data) with time t1 will be interpolated to match time set 2 (t2).
    The returned values, d3, will be data set 1 at time 2.
    No extrapolation will be done; t2's boundaries should encompass those of
    t1.

    Bad data values will be removed prior to interpolation. The `varname`
    kwarg will determine what is considered "bad". Possible varnames include:
    n, t, u[xyz], b[xyz], pos

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

    # Search for bad data values and remove:
    loc = ~np.isfinite(data) | (np.abs(data) >= 1E15)
    if varname in ['n', 'alpha']:
        loc = (loc) | (data > 1E4)
    elif varname == 't':
        loc = (loc) | (data > 1E7)
    elif varname in ['ux', 'uy', 'uz']:
        loc = (loc) | (np.abs(data) > 1E4)
    if args.verbose:
        print(f"Found {loc.sum():05d} bad data points in variable {varname}")

    t1, data = t1[~loc], data[~loc]

    # Remove masked values (if given):
    if type(data) is MaskedArray:
        if type(data.mask) is not bool_:
            d = data[~data.mask]
            t1 = t1[~data.mask]
        else:
            d = data
    else:
        d = data

    # Create interpolator function
    func = interp1d(t1, d, fill_value='extrapolate', **kwargs)

    # Interpolate and return.
    return func(t2)


def read_dscov(fname, fname2=None):
    '''
    Given 1 or 2 of 2 DSCOVR data CDFs (one for the Faraday cup, one for MAG),
    load the data and map to SWMF variables. If only one file name is given,
    the other will be automatically determined based on the first.

    Return a dictionary where each key is an SWMF var name mapped to the
    corresponding value in the DSCOVR data.
    '''
    # Get critical parts of file name:
    result = re.search('dscovr\_h\ds\_(fc|mag)\_(\d+)\_(\d+)(\_cdaweb)?\.cdf',
                       args.file)
    ftype, stime, etime, cdatag = result.groups()

    if args.verbose:
        print(f"Found start/end times in file name of {stime}, {etime}")

    # Set path of files.
    dirname = path.dirname(fname)
    if dirname == '':
        dirname = './'

    # Find partner files; load CDFs.
    if 'fc' in ftype:
        fcp = CDF(fname)
        # Get matching data:
        if not fname2:
            # Build matching name
            basefile = f'dscovr_h?s_mag_{stime}_{etime}{cdatag}.cdf'
            fname2 = glob(path.join(dirname, basefile))[0]
        mag = CDF(fname2)
    elif 'mag' in ftype:
        mag = CDF(fname)
        # Get  matching data:
        if not fname2:
            # Build matching name
            basefile = f'dscovr_h?s_fc_{stime}_{etime}{cdatag}.cdf'
            fname2 = glob(path.join(dirname, basefile))[0]
        fcp = CDF(fname2)
    else:
        raise ValueError('Expected "mag" or "fc" in CDF file name.')

    # Get unified time:
    t_fcp, t_mag = fcp['Epoch'][:], mag['Epoch1'][:]
    time = unify_time(t_fcp, t_mag)

    # Convert coordinates:
    vx, vy, vz = gse_to_gsm(fcp['V_GSE'][:, 0], fcp['V_GSE'][:, 1],
                            fcp['V_GSE'][:, 2], t_fcp)

    # Extract plasma parameters from Faraday cup instrument:
    raw = {'time': time,
           'n': pair(t_fcp, fcp['Np'][:], time, varname='n'),
           't': pair(t_fcp, fcp['THERMAL_TEMP'][:], time, varname='t'),
           'ux': pair(t_fcp, vx, time, varname='ux'),
           'uy': pair(t_fcp, vy, time, varname='uy'),
           'uz': pair(t_fcp, vz, time, varname='uz')}

    # Rotate magnetic field from GSE to GSM:
    bx, by, bz = gse_to_gsm(mag['B1GSE'][:, 0], mag['B1GSE'][:, 1],
                            mag['B1GSE'][:, 2], t_mag)

    # Pair high-time resolution mag data to SWEPAM time:
    for b, bval in zip(['bx', 'by', 'bz'], [bx, by, bz]):
        raw[b] = pair(t_mag, bval, time, varname=b)

    return raw


def read_wind(fname, fname2=None):
    '''
    Given 1 or 2 of 2 Wind data CDFs (one for SWEPAM, one for MAG), load the
    data and map to SWMF variables. If only one file name is given, the
    other will be automatically determined based on the first.

    Return a dictionary where each key is an SWMF var name mapped to the
    corresponding value in the ACE data.
    '''
    # Get critical parts of file name:
    result = re.search('wi(nd)?\_h\ds\_(mfi|swe)\_(\d+)' +
                       '\_(\d+)(\_cdaweb)?\.cdf', fname)

    # If the above doesn't match, use single-file alternative:
    if result is None:
        if args.verbose:
            print("Single-file WIND data detected...")
        # cdf_vars = ['BX', 'BY', 'BZ', 'VX', 'VY', 'VZ', 'Np', 'TEMP']
        # Open solar wind data and load variables into a dictionary:
        obs = CDF(fname)
        raw = {'bx': obs['BX'][...], 'by': obs['BY'][...],
               'bz': obs['BZ'][...], 'ux': obs['VX'][...],
               'uy': obs['VY'][...], 'uz': obs['VZ'][...],
               'n': obs['Np'][...], 't': obs['TEMP'][...],
               'time': obs['Epoch'][...]}
        return raw

    # Otherwise, two-file system. Use names to determine file (MFI vs. SWE)
    has_nd, ftype, stime, etime, cdatag = result.groups()

    if args.verbose:
        print("Two-file WIND data detected.")
        print(f"Found start/end times in file name of {stime}, {etime}")

    # Set path of files.
    dirname = path.dirname(fname)
    if dirname == '':
        dirname = './'

    # Find partner files; load CDFs.
    if 'swe' in ftype:
        swe = CDF(fname)
        # Get matching data:
        if not fname2:
            # Build matching name
            basefile = f'wi*_h?s_mfi_{stime}_{etime}{cdatag}.cdf'
            fname2 = glob(path.join(dirname, basefile))[0]
        mag = CDF(fname2)
    elif 'mfi' in ftype:
        mag = CDF(fname)
        # Get  matching data:
        if not fname2:
            # Build matching name
            basefile = f'wi*_h?s_swe_{stime}_{etime}{cdatag}.cdf'
            fname2 = glob(path.join(dirname, basefile))[0]
        swe = CDF(fname2)
    else:
        raise ValueError('Expected "mfi" or "swe" in CDF file name.')

    # Create unified time:
    t_swe, t_mag = swe['Epoch'][:], mag['Epoch'][:]
    time = unify_time(t_swe, t_mag)

    # Convert coordinates:
    vx, vy, vz = gse_to_gsm(swe['Proton_VX_moment'][:],
                            swe['Proton_VY_moment'][:],
                            swe['Proton_VZ_moment'][:], swe['Epoch'][:])

    # Convert temperature
    temp = (mp/(2*kboltz))*(swe['Proton_W_moment'][:]*1000)**2

    # Extract plasma parameters
    raw = {'time': time,
           'n': pair(t_swe, swe['Proton_Np_moment'][:], time, varname='n'),
           't': pair(t_swe, temp, time, varname='t'),
           'ux': pair(t_swe, vx, time, varname='ux'),
           'uy': pair(t_swe, vy, time, varname='uy'),
           'uz': pair(t_swe, vz, time, varname='uz')}

    # Extract magnetic field parameters:
    for i, b in enumerate(['bx', 'by', 'bz']):
        raw[b] = pair(t_mag, mag['BGSM'][:, i], raw['time'], varname=b)

    # Optional values: Update with more experience w/ wind...
    # if 'alpha_ratio' in swe:
    #     raw['alpha'] = swe['alpha_ratio'][...]
    if 'PGSM' in mag:
        raw['pos'] = pair(t_mag, mag['PGSM'][:, 0],
                          raw['time'], varname='pos') * RE

    return raw


def read_ace(fname, fname2=None):
    '''
    Given 1 of 2 ACE data CDFs (one for SWEPAM, one for MAG), find the matching
    CDF and load the data. Map to SWMF variables.

    Return a dictionary where each key is an SWMF var name mapped to the
    corresponding value in the ACE data.
    '''

    # Get critical parts of file name:
    result = re.search('ac\_h\ds\_(mfi|swe)\_(\d+)\_(\d+)(\_cdaweb)?\.cdf',
                       args.file)
    ftype, stime, etime, cdatag = result.groups()
    if args.verbose:
        print(f"Found start/end times in file name of {stime}, {etime}")

    # Set path of files.
    dirname = path.dirname(fname)
    if dirname == '':
        dirname = './'

    # Find partner files.
    if 'swe' in ftype:
        swe = CDF(fname)
        # Open matching file:
        if not fname2:
            # Build matching name
            basefile = f'ac_h?s_mfi_{stime}_{etime}{cdatag}.cdf'
            fname2 = glob(path.join(dirname, basefile))[0]
        mag = CDF(fname2)
    elif 'mfi' in ftype:
        mag = CDF(fname)
        # Open matching file:
        if not fname2:
            # Build matching name
            basefile = f'ac_h?s_swe_{stime}_{etime}{cdatag}.cdf'
            fname2 = glob(path.join(dirname, basefile))[0]
        swe = CDF(fname2)
    else:
        raise ValueError('Expected "mfi" or "swe" in CDF file name.')

    # Get unified time:
    t_swe, t_mag = swe['Epoch'][:], mag['Epoch'][:]
    time = unify_time(t_swe, t_mag)

    # Extract plasma parameters from SWEPAM
    raw = {'time': time,
           'n': pair(t_swe, swe['Np'][:], time, varname='n'),
           't': pair(t_swe, swe['Tpr'][:], time, varname='t'),
           'ux': pair(t_swe, swe['V_GSM'][:, 0], time, varname='ux'),
           'uy': pair(t_swe, swe['V_GSM'][:, 1], time, varname='uy'),
           'uz': pair(t_swe, swe['V_GSM'][:, 2], time, varname='uz')}

    # Optional values:
    if 'SC_pos_GSM' in mag:
        raw['pos'] = pair(t_mag, mag['SC_pos_GSM'][:, 0], time, 'pos')
    if 'alpha_ratio' in swe:
        raw['alpha'] = pair(t_swe, swe['alpha_ratio'][:], time, 'alpha')

    # Pair high-time resolution mag data to SWEPAM time:
    for i, b in enumerate(['bx', 'by', 'bz']):
        raw[b] = pair(t_mag, mag['BGSM'][:, i], time, varname=b)

    return raw


# ## Begin main script ## #
# Look at filename; determine source and convert data.
if args.verbose:
    print('Reading the following files:')
    print(f'\tFile 1: {args.file}')
    print(f'\tFile 2: {args.file2}')
if (args.file[:2] == 'ac') and ('mfi' in args.file or 'swe' in args.file):
    if args.verbose:
        print('ACE data file detected.')
    raw = read_ace(args.file, args.file2)
elif args.file[:2] == 'wi':
    if args.verbose:
        print('Wind data file detected.')
    raw = read_wind(args.file, args.file2)
elif args.file[:6] == 'dscovr':
    if args.verbose:
        print('DSCOVR data file detected.')
    raw = read_dscov(args.file, args.file2)

# Create seconds-from-start time array:
tsec = np.array([(t - raw['time'][0]).total_seconds() for t in raw['time']])

# Get S/C distance. If not in file, use approximation.
if 'pos' in raw:
    if args.verbose:
        print('S/C location found! Using dynamic location.')
    raw['X'] = raw['pos'] - bound_dist
else:
    if args.verbose:
        print('S/C location NOT found, using static L1 distance.')
    raw['X'] = l1_dist - bound_dist

# Apply velocity smoothing as required
if args.verbose:
    print(f'Applying smoothing using a {args.smoothing} window size.')
velsmooth = medfilt(raw['ux'], args.smoothing)

# Shift time: distance/velocity = timeshift (negative in GSM coords)
if args.tshift >= 0:
    if args.verbose:
        print(f'Using a STATIC timeshift of {args.tshift} minutes.')
        shift = np.zeros(raw['time'].size) + args.tshift * 60.0
else:
    if args.verbose:
        print('Using BALLISTIC propagation.')
    shift = raw['X']/velsmooth  # Time shift per point.
# Apply shift to times.
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
imfout.attrs['header'].append(f'Source data: {args.file}\n')
if args.tshift >= 0:
    imfout.attrs['header'].append('Propagted from L1 to upstream ' +
                                  'BATS-R-US boundary using a constant ' +
                                  f'time delay of {args.tshift} minutes.\n')
else:
    imfout.attrs['header'].append('Ballistically propagted from L1 to ' +
                                  'upstream BATS-R-US boundary\n')
imfout.attrs['header'].append('\n')
imfout.attrs['coor'] = 'GSM'
imfout.attrs['satxyz'] = [np.mean(raw['X'])/RE, 0, 0]
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
