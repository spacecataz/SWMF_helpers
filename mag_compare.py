#!/usr/bin/env python3
'''
This script will create a set of data-model comparisons between a given
SWMF magnetometer output file and any matching magnetometers (by 3 letter IAGA
code) found in a given SuperMag ascii file.  Non matching stations are skipped.

Comparisons will show 3-component data-model comparisons in magnetic
north-east-down coordinates.

Resulting files are placed in a new folder whose name defaults to
"mag_compares".  Plots are saved as "XXX.png" where XXX is the 3-letter IAGA
code.

Requirements
------------
This package requires spacepy version>=0.2 and the SuperMag python reader
found here: https://github.com/spacecataz/supermag
Note that SuperMag tends to change file formats impulsively and often, so
watch for updates to account for this.

Examples
--------
mag_compare.py supermag.txt magnetometers.mag

To-Do:
-------
"Info" flag to print off station info; match in obs file
"--station/-s" option to explicity choose stations to plot.
Set alpha based on number of stations
Set custom time range.
'''

import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Start by configuring the argparser:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("obs", type=str, metavar='observations',
                    help="Path of SuperMag observations file.")
parser.add_argument("mod", nargs='+', metavar='model',
                    help="Path(s) of SWMF model results 'virtual magnetometer'" +
                    " file(s).  Can specify as many as needed.")
parser.add_argument("-l", "--labels", nargs='+', help='Set legend labels for each '+
                    'model included in the comparison')
parser.add_argument("-o", "--outdir", default='mag_compares', help="Set " +
                    "output directory for plots.  Defaults to ./mag_compares")
parser.add_argument("-z", "--horizontal", default=True, help="Instead of dBdown, " +
                    "plot the horizontal component of dB. Requires calculation of " +
                    "the H component in the SuperMAG reader.")
parser.add_argument("--debug", default=False, action='store_true',
                    help="Turn on debugging mode.")
args = parser.parse_args()


# Post-argument imports:
import datetime as dt
import matplotlib.pyplot as plt

from spacepy.pybats import bats
from spacepy.plot import style, applySmartTimeTicks

from supermag import SuperMag, read_statinfo

# Turn on Spacepy plot styles:
style()

# Turn on interactive plotting mode when debug is set:
#if args.debug: plt.ion()

# Load station info:
if args.debug: print('\tLoading magnetometer station info...')
mag_info = read_statinfo()

def comp_mag(name, obs, mod, labels=None, h=False, interactive=False):
    '''
    Given a magnetometer with station name "name", compare the
    data and model together and save the plot.

    Arguments:
    ----------
    name : str
        The 3-letter station code, e.g., "BOU".

    obs : supermag.SuperMag
        A SuperMag-class object containing the observations loaded from a
        SuperMag ascii file.

    mod : list
        A list of spacepy.pybats.MagFile objects to be plotted.

    labels : list
        A list of labels for the plot legend.

    h: bool
        Whether or not to plot the horizontal component instead of the 
        vertical component.
    '''

    # Get number of models to include in this comparison:
    nmod = len(mod)
    
    # Create default set of labels if argument not present:
    if not labels:
        labels = [f'Simulation {i+1}' for i in range(nmod)]

    # Ensure we have at least as many labels as models:
    if len(labels)<nmod:
        for i in range(nmod-len(labels)):
            labels.append(f'Simulation {nmod-i}')
        
    # Set alpha using number of simulations.
    #alpha = 1./nmod -- think of another function.
    alpha = 1.-.3*(nmod>1)
    
    # Get time range that spans all models:
    trange = [dt.datetime(3000,1,1),dt.datetime(1,1,1)]
    for m in mod:
        trange[0]=min(trange[0],m['time'][0])
        trange[1]=max(trange[1],m['time'][-1])
        
    # Create a figure with 3 subplots
    fig = plt.figure(figsize=(10, 10))
    a1, a2, a3 = fig.subplots(3, 1)

    if not h:
        comp1 = 'ned' # model
        comp2 = 'xyz' # obs
    else:
        comp1 = 'neh'
        comp2 = 'xyH'

    # Loop over field component; plot data v. model for each
    for x1, x2, ax in zip(comp1, comp2, (a1, a2, a3)):
        # Plot model & data:
        if name in obs:
            ax.plot(obs['time'], obs[name]['b'+x2], 'k',
                    label='SuperMag Obs.', lw=1.5)
        else:
            ax.plot([],[], 'k', lw=1.5, label='SuperMag Obs.')

        for m,l in zip(mod, labels):
            if name in m:
                if h:
                    # Calculate dBh for the model:
                    m.calc_h()
                
                ax.plot(m[name]['time'], m[name]['dB'+x1], label=l, alpha=alpha)
                ax.set_ylabel(f'$\\Delta B_{x2}$ ($nT$)')

        # Adjust
        applySmartTimeTicks(ax, trange, dolabel=ax == a3)

    # Add legend/title
    a1.legend(loc='best')

    # Build title for our plot using mag_info (if available):
    if name in mag_info:
        info = mag_info[name]
        title = f"{info['station-name']} ({name}) " + \
            f"(mlat={info['aacgmlat']:.1f}$^{{\\circ}}$, " + \
            f" mlon={info['aacgmlon']:.1f}$^{{\\circ}}$)"
    else:
        title = f"Station {name}"
    a1.set_title(title)

    fig.tight_layout()

    if not interactive:
        fig.savefig(f'{args.outdir}/{name}.png')
        plt.close('all')
    else:
        plt.show()


        
# Create output directory
if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)

# Open our data:
if args.debug: print('\tReading observational data file...')
obs = SuperMag(args.obs)

# Check that the observations have the H component:
if args.horizontal:
    if 'bH' not in obs['EYR'].keys():
        print('No bH present in the observations.')
        print('Be sure to turn on calc_H in __init__.py')
        print('in the supermag package (line 336).')
        sys.exit()

if args.debug: print('\tReading model data file(s)...')
mod = [bats.MagFile(x) for x in args.mod]

# Set station list by finding set of all existing
# magnetometers in all model files.
maglist = []
for m in mod:
    for namemag in m.attrs['namemag']:
        if namemag not in maglist: maglist.append(namemag)

# Count total number of magnetometers:
nStats = len(namemag)

# Get list of magnetometers
for i, station in enumerate(maglist):
    print(f'Working on station {i} of {nStats}: {station}')
    if station in obs:
        comp_mag(station, obs, mod, labels=args.labels,
                 h=args.horizontal, interactive=args.debug)
        # If in debug mode, only show 1 plot.
        if args.debug:
            break
    else:
        print('\t...no match.')
