#!/usr/bin/env python
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
"Interactive mode" to show plots on screen/no save.
"--station/-s" option to explicity choose stations to plot.
'''

import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from spacepy.pybats import bats
from spacepy.plot import style, applySmartTimeTicks

from supermag import SuperMag

# Start by configuring the argparser:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("obs", type=str, metavar='observations',
                    help="Path of SuperMag observations file.")
parser.add_argument("mod", type=str, metavar='model',
                    help="Path of SWMF model results 'virtual magnetometer'" +
                    " file.")
parser.add_argument("-o", "--outdir", default='mag_compares', help="Set " +
                    "output directory for plots.  Defaults to ./mag_compares")
args = parser.parse_args()

# Turn on Spacepy plot styles:
style()


def comp_mag(name, obs, mod, interactive=False):
    '''
    Given a magnetometer with station name "name", compare the
    data and model together and save the plot.
    '''
    import matplotlib.pyplot as plt

    # Create a figure with 3 subplots
    fig = plt.figure(figsize=(10, 10))
    a1, a2, a3 = fig.subplots(3, 1)

    # Loop over field component; plot data v. model for each
    for x1, x2, ax in zip('ned', 'xyz', (a1, a2, a3)):
        # Plot model & data:
        ax.plot(mod['time'], mod['dB'+x1], label='SWMF')
        ax.plot(obs['time'], obs[name]['b'+x2], label='SuperMag')
        ax.set_ylabel(f'$\\Delta B_{x2}$ ($nT$)')
        
        # Adjust
        applySmartTimeTicks(ax, mod['time'], dolabel=ax == a3)

    # Add legend/title
    a1.legend(loc='best')
    a1.set_title(f'Magnetometer {name}')

    fig.tight_layout()

    if not interactive:
        fig.savefig(f'{args.outdir}/{name}.png')
        plt.close('all')


# Create output directory
if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)

# Open our data:
obs = SuperMag(args.obs)
mod = bats.MagFile(args.mod)

nStats = len(mod.attrs['namemag'])

# Get list of magnetometers
for i, station in enumerate(mod.attrs['namemag']):
    print(f'Working on station {i} of {nStats}: {station}')
    if station in obs: 
        comp_mag(station, obs, mod[station])
    else:
        print('\t...no match.')
