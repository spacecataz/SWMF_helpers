#!/usr/bin/env python
'''
A script to calculate SMR from SWMF output
and plot it in comparison to the actual
SuperMAG indices. Outputs a pkl file containing
the calculated indexes, as well
as time (datetimes), all in numpy arrays.

Currently only works with mag_grid*.outs 
files, so be sure to concat_outs.py first!

~~~~~~~~~~~~
TO DO: add adjustable time range
'''

import os
import sys
import pickle
from glob import glob
from argparse import ArgumentParser

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

from spacepy.pybats import bats
from spacepy.plot import style, applySmartTimeTicks

import utils

# Turn on spacepy style for plot
style()

parser  = ArgumentParser(description=__doc__)
parser.add_argument("mod", metavar='model', help="Path of SWMF model " +
                    "results 'virtual magnetometer' file.")

args = parser.parse_args()

# Open file and create bats MagGridFile object:
mag = bats.MagGridFile(args.mod)

# Define some strings for the current contributions:
contribs = ['', 'Mhd', 'Fac', 'Hal', 'Ped']
    
def calc_smr(mag, contribs):
    '''
    Takes a Bats mag_grid*.outs file and uses it to 
    calculate the simulated global and regional SMR, 
    split by current contribution. Returns 
    a dictionary of numpy arrays.    
    '''

    # Containers for regional indexes:
    regional = {}
    for i in range(4):
        for c in contribs:
            regional['smr'+str(i*6)+c] = np.zeros(mag.attrs['nframe'])

    times = mag.attrs['times']

    # Container for global indexes:
    globalsmr = {}

    for i in range(mag.attrs['nframe']):
        mag.switch_frame(i)

        # Calculate the MLTs:
        mlts = utils.calc_mlt(mag['LonSm'], mag.attrs['time'])
        # Correct negative MLTs:
        mlts[mlts<0.] += 24.
        
        # Select only grid points between -50 and 50 MLAT:
        latcond = (mag['LatSm'] < 50.) & (mag['LatSm'] > -50.)
        for c in contribs:
            dBn = ma.MaskedArray(mag['dBn'+c], mask=~latcond)
            mlats = ma.MaskedArray(mag['LatSm'], mask=~latcond)

            # Correct for mlat:
            Ncorr = np.divide(dBn, np.cos(np.deg2rad(mlats)))

            # Bin by MLT:
            midcond = (mlts<=3.) | (mlts>21.)
            dawncond = (mlts<=9.) & (mlts>3.)
            nooncond = (mlts<=15.) & (mlts>9.)
            duskcond = (mlts<=21.) & (mlts>15.)

            regional['smr0'+c][i] = np.mean(Ncorr[midcond])
            regional['smr6'+c][i] = np.mean(Ncorr[dawncond])
            regional['smr12'+c][i] = np.mean(Ncorr[nooncond])
            regional['smr18'+c][i] = np.mean(Ncorr[duskcond])

    # Calculate global SMR:
    for c in contribs:
        globalsmr['smr'+c] = (regional['smr0'+c]+regional['smr6'+c] \
                              +regional['smr12'+c]+regional['smr18'+c]) / 4
        
    return times, regional, globalsmr

# Calculate the grid-based indexes:
times, regional, globalsmr = calc_smr(mag, contribs)
result = {'times': times, 'regional': regional, 'global': globalsmr}

with open('grid_smr.pkl', 'wb') as f:
    pickle.dump(result, f)
    print('Data saved as ./grid_smr.pkl')

# Plot the calculations as a check:
fig, axs = plt.subplots(5, 1, figsize=(14,12), sharex=True)
for j, c in enumerate(contribs):
    # Plot the globals:
    axs[j].plot(times, globalsmr['smr'+c], color='k', label='global')
    axs[j].set_ylabel('SMR '+c+' (nT)')

    # Plot the regionals:
    for k, s in enumerate(['smr'+str(i*6)+c for i in range(4)]):
        axs[j].plot(times, regional[s], label='smr'+str(k*6))


fig.legend(['SMR', 'SMR00', 'SMR06', 'SMR12', 'SMR18'],
           loc='outside right')
fig.suptitle('Low Latitude Indexes Calculated From mag_grid*.outs')

applySmartTimeTicks(axs[-1], times, dolabel=True)

plt.savefig('grid_smr.png')
plt.show()

