#!/usr/bin/env python3
'''
# Solar Wind Toolbox

This module supports solar-wind related tools within this repository.
It contains commonly used functions for fetching, cleaning, and preparing data
for use within the SWMF.

## Data Fetching/Loading Utilities


'''

from glob import glob
from os import path
import re
from datetime import datetime, timedelta

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import medfilt
from matplotlib.lines import Line2D

from hapiclient import hapi
from spacepy.pycdf import CDF
from spacepy.datamodel import dmarray
from spacepy.pybats import ImfInput

# Information for HAPI fetching of data:
hapiserv = 'https://cdaweb.gsfc.nasa.gov/hapi'

# Declare important constants
RE = 6371              # Earth radius in kmeters.
l1_dist = 1495980      # L1 distance in km.
kboltz = 1.380649E-23  # Boltzmann constant, J/K
mp = 1.67262192E-27    # Proton mass in Kg
bound_dist = 32 * RE   # Distance to BATS-R-US upstream boundary from Earth.

# Set var names and units used by the SWMF.
swmf_vars = ['bx', 'by', 'bz', 'ux', 'uy', 'uz', 'n', 't']
units = {v: u for v, u in zip(swmf_vars, 3*['nT']+3*['km/s']+['cm-3', 'K'])}


def convert_hapi_t(time):
    '''
    Convert HAPI time stamps into arrays of datetime objects.
    '''

    from dateutil.parser import parse

    return np.array([parse(t.decode('utf-8')) for t in time])


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


def pair(time1, data, time2, varname=None, verbose=False, **kwargs):
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
    if verbose:
        print(f"Found {loc.sum():05d} bad data points in variable {varname}")

    t1, data = t1[~loc], data[~loc]
    if data.size == 0:
        raise ValueError('No valid data points in this time range.')

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


def read_ace_cdf(fname, fname2=None, outname=None, verbose=False):
    '''
    Given 1 of 2 ACE data CDFs (one for SWEPAM, one for MAG), find the matching
    CDF and load the data. Build an SWMF input file object and return.

    ACE spacecraft CDFs come in pairs: a SWEPAM file with plasma information,
    spacecraft position, etc. and an MFI file with magnetic field data. If you
    specify either file as the input argument, this script will attempt to find
    the matching file in the same directory. On CDAWeb, select `AC_H0_MFI` and
    `AC_H0_SWE` datasets and obtain values in GSM coordinates.
    ACE spacecraft CDFs must contain the following variables:
    | Variable Name(s) | Description                                          |
    |------------------|------------------------------------------------------|
    | SC_pos_GSM       | Distance from Earth (km) in GSM/GSE (optional)       |
    | alpha_ratio      | Alpha to proton ratio (not currently used)           |
    | BGSM             | Vector magnetic field (nT) data in GSM coordinates   |
    | V_GSM            | Vector solar wind velocity (km/s) in GSM coordinates |
    | Np               | Proton number density (1/ccm)                        |
    | Tpr              | Plasma temperature (K)                               |
    | Epoch            | CDF-formatted universal time.                        |

    Parameters
    ----------
    fname : str
        The file path/name for 1 of 2 CDF files for creating an SWMF input
        file based on downloaded ACE CDFs. Can either be a SWEPAM or MAG
        data set.
    fname2 : str, defaults to None
        The file path/name for the second CDF file. Specify this argument if
        the function is not automatically finding the file that pairs with
        the first file (`fname`)
    outname : str, defaults to None
        If given, sets the output file name in the resulting SWMF object and
        saves the file to disk.
    verbose : bool, defaults to False
        Activate verbose mode.

    Returns
    -------
    swout : spacepy.pybats.ImfInput object
        The data converted into an SWMF-formatted IMF object. Two additional
        keys will be added if found in the CDFs: 'alpha' (the alpha ratio) and
        'pos' (the spacecraft distance from the Earth in km).

    '''

    # Get critical parts of file name:
    result = re.search('ac\_h\ds\_(mfi|swe)\_(\d+)\_(\d+)(\_cdaweb)?\.cdf',
                       fname)
    ftype, stime, etime, cdatag = result.groups()
    if verbose:
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

    # Create IMF object:
    npts = raw['time'].size
    swout = ImfInput(filename=outname, load=False, npoints=npts)

    for v in swmf_vars:
        swout[v] = dmarray(raw[v], {'units': units[v]})
    swout['time'] = raw['time']
    swout['alpha'] = dmarray(raw['alpha'], {'units': 'ratio'})
    swout['pos'] = dmarray(raw['pos'], {'units': 'km'})

    # Build header:
    avgdist = raw['pos'].mean() / RE
    swout.attrs['header'].append(f'Source data: \n\t{fname}\n\t{fname2}\n')
    swout.attrs['header'].append(f'File created on {datetime.now()}\n')
    swout.attrs['header'].append(f'SC mean distance from Earth: {avgdist}RE\n')
    swout.attrs['header'].append('\n')
    swout.attrs['coor'] = 'GSM'

    # Save data as necessary:
    if outname:
        swout.write()

    return swout


def fetch_ace_hapi(tstart, tend, outname=None, verbose=False):
    '''
    Fetch ACE solar wind data from the CDAWeb HAPI server and convert to an
    SWMF ImfInput object.

    Bad data flags are removed and coordinates rotated to GSM, but no other
    processing (including propagation from L1) is performed

    Parameters
    ----------
    tstart, tend : datetime.datetime objects
        The start and end time of the interval to fetch.
    outname : str, defaults to None
        If given, sets the output file name in the resulting SWMF object and
        saves the file to disk.
    verbose : bool, defaults to False
        Activate verbose mode.

    Returns
    -------
    swout: spacepy.pybats.ImfInput
        The resulting ImfInput object.
    '''

    swedat, magdat = 'AC_H0_SWE', 'AC_H0_MFI'
    swevar = 'Np,Tpr,alpha_ratio,V_GSM'
    magvar = 'BGSM,SC_pos_GSM'

    t1, t2 = tstart.isoformat(), tend.isoformat()
    swe, meta = hapi(hapiserv, swedat, swevar, t1, t2)
    mag, meta = hapi(hapiserv, magdat, magvar, t1, t2)

    if swe.size <= 1 or mag.size <= 1:
        raise ValueError('No data for SWE or MFI')

    # Get unified time:
    t_swe, t_mag = convert_hapi_t(swe['Time']), convert_hapi_t(mag['Time'])
    time = unify_time(t_swe, t_mag)

    # Extract plasma parameters from SWEPAM
    raw = {'time': time,
           'n': pair(t_swe, swe['Np'][:], time, varname='n'),
           't': pair(t_swe, swe['Tpr'][:], time, varname='t'),
           'ux': pair(t_swe, swe['V_GSM'][:, 0], time, varname='ux'),
           'uy': pair(t_swe, swe['V_GSM'][:, 1], time, varname='uy'),
           'uz': pair(t_swe, swe['V_GSM'][:, 2], time, varname='uz')}

    # Optional values:
    if 'SC_pos_GSM' in mag.dtype.names:
        raw['pos'] = pair(t_mag, mag['SC_pos_GSM'][:, 0], time, 'pos')
    if 'alpha_ratio' in swe.dtype.names:
        raw['alpha'] = pair(t_swe, swe['alpha_ratio'][:], time, 'alpha')

    # Pair high-time resolution mag data to SWEPAM time:
    for i, b in enumerate(['bx', 'by', 'bz']):
        raw[b] = pair(t_mag, mag['BGSM'][:, i], time, varname=b)

    # Create IMF object:
    npts = raw['time'].size
    swout = ImfInput(filename=outname, load=False, npoints=npts)

    for v in swmf_vars:
        swout[v] = dmarray(raw[v], {'units': units[v]})
    swout['time'] = raw['time']
    swout['alpha'] = dmarray(raw['alpha'], {'units': 'ratio'})
    swout['pos'] = dmarray(raw['pos'], {'units': 'km'})

    # Build header:
    avgdist = raw['pos'].mean() / RE
    swout.attrs['header'].append('Data obtained from CDAWeb HAPI\n')
    swout.attrs['header'].append(f'File created on {datetime.now()}\n')
    swout.attrs['header'].append(f'SC mean distance from Earth: {avgdist}RE\n')
    swout.attrs['header'].append('\n')
    swout.attrs['coor'] = 'GSM'

    # Save data as necessary:
    if outname:
        swout.write()

    return swout


def fetch_wind_hapi(tstart, tend, outname=None, verbose=False):
    # Variables for interfacing with CDAweb's HAPI server:

    swedat, magdat = 'WI_H1S_SWE', 'WI_H0_MFI@1'
    swevar = ''
    magvar = 'B3GSM'

    swe, meta = hapi(hapiserv, swedat, swevar, tstart, tend)
    mag, meta = hapi(hapiserv, magdat, magvar, tstart, tend)


def fetch_dscovr_hapi(tstart, tend, outname=None, verbose=False):
    pass


def fetch_hapi(source, tstart, tend, outname=None, verbose=False):
    '''
    Convenience function for fetching solar wind data from CDAWeb via HAPI.

    Parameters
    ----------
    source : str
        Select source: 'ace', 'dscovr', 'wind'
    tstart, tend : datetime.datetime objects
        The start and end time of the interval to fetch.
    outname : str, defaults to None
        If given, sets the output file name in the resulting SWMF object and
        saves the file to disk.
    verbose : bool, defaults to False
        Activate verbose mode.

    Returns
    -------
    swout: spacepy.pybats.ImfInput
        The resulting ImfInput object.
    '''

    funcs = {
        'ace': fetch_ace_hapi,
        'dscovr': fetch_dscovr_hapi,
        'wind': fetch_wind_hapi
             }

    return funcs[source](tstart, tend, outname, verbose)


def l1_propagate(swfile, outfile=None, shift=-1.0, smoothwin=1,
                 verbose=False):
    '''
    Given an SWMF ImfInput object or file with data assumed to be at L1,
    propagate to the SWMF upstream boundary position.

    The default behavior is to propagate ballistically, removing points that
    are overtaken by subsequent points. If a time shift is given, a constant
    shift is used instead.

    The position of the spacecraft defaults to 1.49Mkm upstream. However, if
    the key 'pos' is found in the input data and it contains the SC distance
    from the Earth in km, a dynamic distance is used.

    For noisy data, it may be desireable to smooth the solar wind velocity.
    The argument `smoothwin` activates median smoothing of velocity and can be
    set to the window size (in number of points) for the median filter.

    Parameters
    ----------
    swfile : str or spacepy.pybats.ImfInput object
        The data to propagate, either as an input filename or an already
        opened data object.
    outfile : str, defaults to None
        The name of the output file. If given, data is saved to disk after
        propagtion.
    shift : float, defaults to -1
        The time shift given in minutes. If -1, use dynamic ballistic
        propagation (default behavior).
    smoothwin : int, defaults to 1
        Sets the window size for median filtering of the solar wind velocity.
    verbose : bool, defaults to False
        Activate verbose output.
    '''

    # Get data.
    if type(swfile) is ImfInput:
        raw = swfile
    else:
        raw = ImfInput(swfile)

    # Create seconds-from-start time array:
    tstart = raw['time'][0]
    tsec = np.array([(t - tstart).total_seconds() for t in raw['time']])

    # Get S/C distance. If not in file, use approximation.
    if 'pos' in raw:
        if verbose:
            print('S/C location found! Using dynamic location.')
        raw['X'] = raw['pos'] - bound_dist
    else:
        if verbose:
            print('S/C location NOT found, using static L1 distance.')
        raw['X'] = l1_dist - bound_dist

    # Apply velocity smoothing as required
    if verbose:
        print(f'Applying smoothing using a {smoothwin} window size.')
    velsmooth = medfilt(raw['ux'], smoothwin)

    # Shift time: distance/velocity = timeshift (negative in GSM coords)
    if shift >= 0:
        if verbose:
            print(f'Using a STATIC timeshift of {shift} minutes.')
            shift = np.zeros(raw['time'].size) + shift * 60.0
    else:
        if verbose:
            print('Using BALLISTIC propagation.')
        lags = raw['X']/velsmooth  # Time shift per point.
    # Apply shift to times.
    tshift = np.array([t1 - timedelta(seconds=t2) for t1, t2 in
                       zip(raw['time'], lags)])

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
    if verbose:
        print(f'Removing "overtaken" points {len(discard)}/{tsec.size} total.')

    # Create new IMF object and populate with propagated values.
    # Use the information above to throw out overtaken points.
    if outfile is not None:
        outfile = outfile + '.dat'*(outfile[-4:] != '.dat')
    imfout = ImfInput(outfile, load=False, npoints=len(keep))
    for v in swmf_vars:
        imfout[v] = dmarray(raw[v][keep], {'units': units[v]})
    imfout['time'] = tshift[keep]
    imfout.attrs['header'] = raw.attrs['header']

    if shift >= 0:
        imfout.attrs['header'].append('Propagted from L1 to upstream ' +
                                      'BATS-R-US boundary using a constant ' +
                                      f'time delay of {shift} minutes.\n')
    else:
        imfout.attrs['header'].append('Ballistically propagted from L1 to ' +
                                      'upstream BATS-R-US boundary\n')
    imfout.attrs['header'].append('\n')
    imfout.attrs['coor'] = 'GSM'

    # Write to disk:
    if outfile:
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

    figname = outfile[:-4] if outfile else raw.attrs['file'][:-4]
    fig.savefig(figname + '_prop_info.png')
