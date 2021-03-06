#!/usr/bin/env python3

'''
For a SWMF simulation for space weather at Earth applications, make a 
series of ionospheric summary plots.

To use this script, call it from the base directory of a completed and
post-processed SWMF results directory.

'''

import os
from sys import stdout
from dateutil.parser import parse
from argparse import ArgumentParser

# Initialize argument parser:
parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--imffile", default='',
                    help="Path to relevant SWMF-formatted IMF input file.  "+
                    "Default action is to search in current directory for "+
                    "imf*dat")
parser.add_argument("-colat", "--colat", type=int, default=40,
                    help="Set co-latitude plot limit for ionosphere " +
                    "contour plots.")

# Handle arguments:
args = parser.parse_args()

from glob import glob

import numpy as np
import matplotlib.pyplot as plt

from spacepy.plot import style, applySmartTimeTicks
from spacepy.pybats import rim, ImfInput

# Switch to spacepy's style:
style()

# Find the IMF file:
imfpath = glob('./imf*.dat')[0] if not args.imffile else args.imffile
imf = ImfInput(imfpath)

# Some constants:
bzcolor = '#3333CC'
bycolor = '#ff9900'
pdcolor = '#CC3300'

outdir='iono_figs/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

def create_figure(iefile):
    # Open IE file, calculate azimuthal current:
    ie = rim.Iono(iefile)
    ie.calc_j()
    ie['n_jphi']/=1000.
    ie['s_jphi']/=1000.
    
    # Create figure:
    fig = plt.figure(figsize=(10,10))
    fig.subplots_adjust(left=.08, bottom=.052, top=.95,
                        right=.962, wspace=.217, hspace=.23)
    # Set common keyword arguments:
    kwargs = {'max_colat':args.colat, 'target':fig, 'add_cbar':True}
    
    # Add IE plots:
    # FACs:
    with plt.style.context('default'):
        zmax = max(np.abs(ie['n_jr']).max(), np.abs(ie['s_jr']).max())
        ie.add_cont('n_jr', loc=431, **kwargs)
        ie.add_cont('s_jr', loc=434, label='', **kwargs)
        # J_phi:
        zmax = max(np.abs(ie['n_jphi']).max(), np.abs(ie['s_jphi']).max())
        ie.add_cont('n_jphi', loc=432, **kwargs)
        ie.add_cont('s_jphi', loc=435, label='', **kwargs)
        # J_phi:
        zmax = max(np.abs(ie['n_phi']).max(), np.abs(ie['s_phi']).max())
        ie.add_cont('n_phi', loc=433, **kwargs)
        ie.add_cont('s_phi', loc=436, label='', **kwargs)
    
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
    applySmartTimeTicks(a5, imf['time'], dolabel=True)

    a4.vlines(ie.attrs['time'], a4.get_ylim()[0], a4.get_ylim()[1], colors='k',
              linestyles='dashed', linewidths=2)
    lim = a5.get_ylim()
    a5.vlines(ie.attrs['time'], lim[0], lim[1], colors='k',
              linestyles='dashed', linewidths=2)
    
    fig.savefig(outdir+f'iono_t{ie.attrs["time"]:%Y%m%d_%H%M%S}.png')
    
    return fig
    
if __name__ == '__main__':
    iefiles = glob('IE/it*.idl')
    if not iefiles: iefiles=glob('IE/it*.idl.gz')

    for ie in iefiles:
        print(f'Working on file {ie}...')
        fig=create_figure(ie)
        plt.close('all')
