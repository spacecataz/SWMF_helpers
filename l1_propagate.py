#!/usr/bin/env python3
'''
Ballistically propagate solar wind parameters from L1 to +32RE upstream for
use in the Space Weather Modeling Framework.

Given a CDF file that contains the requisite values (see below), the time
array is delayed by the time each solar wind "parcel" (point in the file)
would take to travel from the L1 position to the upstream BATS-R-US boundary.
The travel time is taken using the parcel's current Earthward flow speed
(Vx), yielding a dynamic time adjustment.

CURRENT ASSUMPTIONS (to change):
- GSM coordinates for B, V vectors.
- The position of the S/C is not used, a constant L1 distance is employed.
- Works only with Wind CDF data files.
- IMF and solar wind values have same time cadence.

Wind spacecraft CDFs must contain the following variables:
| Variable Name(s) | Description                                            |
|------------------|--------------------------------------------------------|
| BX, BY, BZ       | Vector magnetic field (nT) data in GSM coordinates.    |
| VX, VY, VZ       | Vector solar wind velocity (km/s) in GSM coordinates.  |
| Np               | Proton number density (1/ccm)                          |
| TEMP             | Plasma temperature (K)                                 |
| Epoch            | CDF-formatted universal time.                          |
'''

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
                    help='A Wind dataset in CDF format that contains the ' +
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
args = parser.parse_args()

# Declare important constants
RE = 6371  # Earth radius in kmeters.
l1_dist = 1495980  # L1 distance in km.
bound_dist = 32 * RE  # Distance to BATS-R-US upstream boundary from Earth.
swmf_vars = ['bx', 'by', 'bz', 'ux', 'uy', 'uz', 'n', 't']
cdf_vars = ['BX', 'BY', 'BZ', 'VX', 'VY', 'VZ', 'Np', 'TEMP']
units = {v: u for v, u in zip(swmf_vars, 3*['nT']+3*['km/s']+['cm-3', 'K'])}

# Open solar wind data and load variables into a dictionary:
obs = CDF(args.file)
raw = {'bx': obs['BX'][...], 'by': obs['BY'][...], 'bz': obs['BZ'][...],
       'ux': obs['VX'][...], 'uy': obs['VY'][...], 'uz': obs['VZ'][...],
       'n': obs['Np'][...], 't': obs['TEMP'][...]}
time = obs['Epoch'][...]

# Create seconds-from-start time array:
tsec = np.array([(t - time[0]).total_seconds() for t in time])

# Interpolate over bad data:
if args.verbose:
    print('Removing bad data:')
for v in swmf_vars:
    # Find bad values:
    loc = ~np.isfinite(raw[v])
    if args.verbose:
        print(f'\t{v} has {loc.sum()} bad values.')
    if loc.sum() == 0:
        continue
    interp = interp1d(tsec[~loc], raw[v][~loc])
    raw[v][loc] = interp(tsec[loc])

# ## TEMP ONLY ## #
# loc = time < datetime(2024, 5, 10, 15, 0, 0)
# raw['n'][loc] = medfilt(raw['n'][loc], 31)
# for x in 'xyz':
#     raw['u'+x][loc] = medfilt(raw['u'+x][loc], 31)
# ############### #

# Apply velocity smoothing as required
if args.verbose:
    print(f'Applying smoothing using a {args.smoothing} window size.')
velsmooth = medfilt(raw['ux'], args.smoothing)

# Shift time: distance/velocity = timeshift (negative in GSM coords)
shift = (l1_dist - bound_dist)/velsmooth  # Time shift per point.
tshift = np.array([t1 - timedelta(seconds=t2) for t1, t2 in zip(time, shift)])

# Ensure that any points that are "overtaken" (i.e., slow wind overcome by
# fast wind) are removed. First, locate those points:
keep = [0]
discard = []
lasttime = tshift[0]
for i in range(1, time.size):
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
plotvars = ['BY', 'BZ', 'Np', 'TEMP', 'VX']
for ax, v in zip(fig.axes, plotvars):
    c = ax.get_lines()[0].get_color()
    ax.plot(obs['Epoch'], obs[v], '--', c=c, alpha=.5)
    ax.plot(obs['Epoch'][...][discard], obs[v][...][discard],
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
