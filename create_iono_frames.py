#!/usr/bin/env python3

'''
For a SWMF simulation for space weather at Earth applications, make a
series of ionospheric summary plots that include polar plots of both
northern and southern hemisphere FACs, horizontal currents, and electric
potential over a time series of solar wind values.

By default, the entire simulation time period will be plotted. This can be
changed using the "start" and "end" arguments, which not only changes which
files will be opened and plotted, but also the range of the timeseries
sub-figure.

To use this script, call it from the base directory of a completed and
post-processed SWMF results directory. The following file types should
be present:
- The IMF input file, either specified using the "imffile" argument or
  present working directory and named "imf*.dat".
- "it*.idl" or "it*.idl.gz" files inside of the IE directory.
'''

import os
from glob import glob
from argparse import ArgumentParser

from dateutil.parser import parse
import matplotlib
import matplotlib.pyplot as plt

from spacepy.plot import style, applySmartTimeTicks
from spacepy.pybats import rim, ImfInput

# Initialize argument parser:
parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--imffile", default='',
                    help="Path to relevant SWMF-formatted IMF input file.  " +
                    "Default action is to search in current directory for " +
                    "imf*dat")
parser.add_argument("-v", "--vars", default=['jr', 'jphi', 'phi'], nargs=3,
                    help="Set the values to plot from jr, jphi, phi, sigmap," +
                    " and sigmah.  Must provide 3 values. Default is " +
                    "jr, jphi and phi.")
parser.add_argument("-j", "--maxj", type=float, default=None,
                    help="Set the horizontal current color bar max. " +
                    "Defaults to file's maximum value for each frame.")
parser.add_argument("-f", "--maxfac", type=float, default=None,
                    help="Set the FAC color bar max. " +
                    "Defaults to file's maximum value for each frame.")
parser.add_argument("-p", "--maxpot", type=float, default=None,
                    help="Set electric potential color bar max. " +
                    "Defaults to file's maximum value for each frame.")
parser.add_argument("-s", "--maxcond", type=float, default=None,
                    help="Set conductance color bar max. " +
                    "Defaults to file's maximum value for each frame.")
parser.add_argument("-colat", "--colat", type=int, default=40,
                    help="Set co-latitude plot limit for ionosphere " +
                    "contour plots.")
parser.add_argument('-start', default=None, help="Set the start date and " +
                    "time of the simulation using the format " +
                    "YYYY-MM-DDTHH:MN:SS.")
parser.add_argument('-d', '--debug', action='store_true',
                    help="Turn on debug output.")
parser.add_argument('-cont', '--cont', action='store_true',
                    help="If set, it will 'continue' a previous call by " +
                    "skipping all files where an output PNG already exists." +
                    " This is useful for continuing an aborted session.")
parser.add_argument('-end', default=None,
                    help="Same as -start but for end time")

# Handle arguments:
args = parser.parse_args()
if not args.debug:
    matplotlib.use('Agg')

# Switch to spacepy's style:
style()

# Find the IMF file:
imfpath = glob('./[iI][mM][fF]*.dat')[0] if not args.imffile else args.imffile
imf = ImfInput(imfpath)

# Some plot constants:
bzcolor = '#3333CC'
bycolor = '#ff9900'
pdcolor = '#CC3300'
outdir = 'iono_figs/'

# Set some labels:
labs = {'jphi': r'$J_{\phi}$', 'phi': None, 'sigmah': None, 'sigmap': None,
        'jr': None}
mz = {'jphi': args.maxj, 'sigmah': args.maxcond, 'phi': args.maxpot,
      'jr': args.maxfac, 'sigmap': args.maxcond}

if not os.path.exists(outdir):
    os.mkdir(outdir)


def parse_ie_time(filename):
    '''Get a datetime from an IE file name'''
    from datetime import datetime

    zipped = '.gz' in filename
    return datetime.strptime(filename[-21-zipped*3:-8-zipped*3],
                             '%y%m%d_%H%M%S')


def absmax(ie, val):
    '''Get the maximum absolute value for both hemispheres.'''
    from numpy import abs

    return max(abs(ie['n_'+val]).max(), abs(ie['s_'+val]).max())


def create_figure(iefile, maxpot=None, trng=None, outname='iono_fig.png'):
    try:
        ie = rim.Iono(iefile)
    except (ValueError, IndexError) as error:
        print(f'ERROR opening file {iefile} ({error})')
        return None

    # Calculate azimuthal current:
    if 'n_jx' in ie:
        ie.calc_j()
        if 'n_jphi' not in ie:
            raise ValueError('Horizontal currents not in IE output.')
        ie['n_jphi'] /= 1000.
        ie['s_jphi'] /= 1000.

    # Set colorbar ranges if not already set:
    for x in args.vars:
        if mz[x] is None:
            mz[x] = absmax(ie, x)

    # Create figure:
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(left=.08, bottom=.052, top=.95,
                        right=.962, wspace=.217, hspace=.23)
    # Set common keyword arguments:
    kwargs = {'max_colat': args.colat, 'target': fig, 'add_cbar': True,
              'extend': 'both'}

    # Add IE plots:
    with plt.style.context('default'):
        for i, x in enumerate(args.vars):
            ie.add_cont('n_'+x, loc=431+i, maxz=mz[x], label=labs[x], **kwargs)
            ie.add_cont('s_'+x, loc=434+i, maxz=mz[x], label=labs[x], **kwargs)

    a4 = fig.add_subplot(413)
    a5 = fig.add_subplot(414, sharex=a4)
    # Add IMF data:
    a4.plot(imf['time'], imf['bz'], color=bzcolor, lw=1.25)
    a4.plot(imf['time'], imf['by'], color=bycolor, lw=1.25)
    a4.hlines(0, imf['time'][0], imf['time'][-1], color='k',
              linestyles='dashed')
    a4.legend(['B$_Z$', 'B$_Y$'], loc='upper left')
    a4.set_ylabel('IMF ($nT$)')

    # Add Pdyn:
    a5.plot(imf['time'], imf['pram'], color=pdcolor)
    a5.grid(False, axis='y')
    a5.set_ylabel('P$_{dyn}$ ($nPa$)')
    applySmartTimeTicks(a5, trng, dolabel=True)

    # Turn on vertical lines for current time:
    a4.vlines(ie.attrs['time'], a4.get_ylim()[0], a4.get_ylim()[1], colors='k',
              linestyles='dashed', linewidths=2)
    lim = a5.get_ylim()
    a5.vlines(ie.attrs['time'], lim[0], lim[1], colors='k',
              linestyles='dashed', linewidths=2)

    fig.savefig(outname)

    return fig


if __name__ == '__main__':
    # Get list of files to work on:
    iefiles = glob('IE/it*.idl')
    if not iefiles:
        iefiles = glob('IE/it*.idl.gz')
    iefiles.sort()

    if args.debug:
        print(f"Found {len(iefiles)} files.")

    # Set time limits for plot:
    if args.start:
        # Convert to datetime.
        start = parse(args.start)
    else:
        # Get time from first file:
        if args.debug:
            print(f"Obtaining start time from {iefiles[0]}...")
        start = parse_ie_time(iefiles[0])

    # And again for end time:
    if args.end:
        # Convert to datetime.
        end = parse(args.end)
    else:
        # Get time from first file:
        if args.debug:
            print(f"Obtaining start time from {iefiles[0]}...")
        end = parse_ie_time(iefiles[-1])

    print(f"Processing files from {start} to {end}")

    for ie in iefiles:
        # Check time against limits:
        tnow = parse_ie_time(ie)
        if (tnow < start) or (tnow > end):
            continue

        # Create figure file name:
        outname = outdir+f'iono_t{tnow:%Y%m%d_%H%M%S}.png'

        # Skip if file exits in "continue" mode.
        if args.cont and os.path.exists(outname):
            continue

        # Create figure for current file:
        print(f'Working on file {ie}...')
        fig = create_figure(ie, trng=[start, end], outname=outname)
        if args.debug:
            plt.show()
            break
        else:
            plt.close('all')
