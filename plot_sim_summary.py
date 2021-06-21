#!/usr/bin/env python3

'''
For a SWMF simulation for space weather at Earth applications, make a 
series of summary plots for quick-look analysis.

To use this script, call it from the base directory of a completed and
post-processed SWMF results directory.

'''

import os
from sys import stdout
from dateutil.parser import parse
from argparse import ArgumentParser, RawDescriptionHelpFormatter
# Required defaults:
start = False
end   = False

# Set up argruments:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-m", "--mags", action="store_true",
                    help="Create summary plots of all virtual magnetometers.")
parser.add_argument("--supermag", default='',
                    help="Plot magnetometer observations over model results "+
                    "using values in SuperMag output file given by argument "+
                    "value (e.g, --supermag 'supermag_data_file.dat'.")
parser.add_argument("-i", "--imffile", default='',
                    help="Path to relevant SWMF-formatted IMF input file.  "+
                    "Default action is to search in current directory for imf*dat")
parser.add_argument("-c", "--currents", default='',
                    help="Path to pickle containing integrated currents from IE "+
                    "as calculated by calc_rim_I.py.  Default action is to grab "+
                    "first found pickle file (*.pkl)")
parser.add_argument("-obs", "--obs", help="Fetch and plot observed indices on "+
                    " top of modeled values", action="store_true")
parser.add_argument("--debug", help="Only the first plot is created.  "+
                    "It is shown to the user and not saved.  Additional "+
                    "information is printed to screen.", action="store_true")

# Handle arguments:
args = parser.parse_args()

from glob import glob

import numpy as np
import matplotlib.pyplot as plt

from spacepy.plot import style, applySmartTimeTicks
from spacepy.pybats import bats, rim, ImfInput

style()

##################### PLOT SUBROUTINES ########################
def read_supermag(filename):

    '''
    Read a complicated supermag file and return a dictionary of 'time'
    and numpy arrays of magetometer delta-Bs for each mag in the file.
    '''

    import datetime as dt
    from matplotlib.dates import date2num
    from scipy.interpolate import interp1d
    
    f = open(filename, 'r')

    # Skip header:
    line = f.readline()
    while 'Selected parameters' not in line: line = f.readline()
    head = f.readline()
    
    # Get station list:
    stats  = head.split()[-1].split(',')
    nStats = len(stats)

    # Now, slurp rest of lines and count number of records.
    f.readline() # Skip last header line.
    lines = f.readlines()

    # Get number of lines:
    nTime = int(len(lines)/(nStats+1))

    # Create container of data:
    data = {}
    data['time'] = np.zeros(nTime, dtype=object)
    for s in stats: # Initialize with Bad Data Flag
        data[s] = np.zeros( [3,nTime] )+999999.

    # Read data:
    for j in range(nTime):
        # Get time:
        data['time'][j] = dt.datetime.strptime(
            ''.join(lines.pop(0).split()[:-1]), '%Y%m%d%H%M%S')
        
        # Get values:
        for i in range(nStats):
            if lines[0][:3] not in stats: continue
            parts = lines.pop(0).split()
            data[parts[0]][:,j] = parts[1:4]

    # Filter bad data.
    t = date2num(data['time'])
    for s in stats:
        data[s][data[s]>=999999.] = np.nan

        # Interpolate over bad data:
        for idir in range(3):
            bad = np.isnan(data[s][idir,:])
            good= np.logical_not(bad)
            data[s][idir,bad] = interp1d(t[good], data[s][idir,good],
                                         fill_value='extrapolate')(t[bad])
            
    # Get time in seconds:
    dt = np.array([x.total_seconds() for x in np.diff(data['time'])])
    
    # Calc H and time derivatives:
    for s in stats:
        # Horizontal field magnitude:
        data[s+'_H'] = np.sqrt(data[s][0,:]**2 + data[s][1,:]**2)

        # Get dB_n/dt and dB_e/dt:
        dbn, dbe = np.zeros(dt.size+1),np.zeros(dt.size+1)

        # Central diff:
        dbn[1:-1] = (data[s][0,2:]-data[s][0,:-2])/(dt[1:]+dt[:-1])
        dbe[1:-1] = (data[s][1,2:]-data[s][1,:-2])/(dt[1:]+dt[:-1])
        # Forward diff:
        dbn[0]=(-data[s][0,2]+4*data[s][0,1]-3*data[s][0,0])/(dt[1]+dt[0])
        dbe[0]=(-data[s][1,2]+4*data[s][1,1]-3*data[s][1,0])/(dt[1]+dt[0])
        # Backward diff:
        dbn[-1]=(3*data[s][0,-1]-4*data[s][0,-2]+data[s][0,-3])/(dt[-1]+dt[-2])
        dbe[-1]=(3*data[s][1,-1]-4*data[s][1,-2]+data[s][1,-3])/(dt[-1]+dt[-2])

        # Create |dB/dt|_h:
        data[s+'_dH'] = np.sqrt(dbn**2 + dbe**2)


    return data


def plot_mags():
    '''
    Plot all magnetometers and save to file.
    '''

    path = './mag_summaries/'
    
    # Create output directory as needed:
    if not os.path.exists(path):
        os.mkdir(path)

    # Open magnetometer files:
    mag_path = glob('GM/magnetometers*.mag')
    if not mag_path:
        raise ValueError('Could not find magnetometers in run directory.')
    mags = bats.MagFile(mag_path[0])
    mags.calc_h()
    mags.calc_dbdt()

    # Open observations as necessary:
    do_obs = bool(args.supermag)
    if do_obs:
        data = read_supermag(args.supermag)

    # Loop through magnetometers:
    stations = mags.attrs['namemag']

    for s in stations:
        # shortcut!
        mag = mags[s]
        
        # Create figure & axes:
        fig=plt.figure(figsize=(10,7.5))
        fig.subplots_adjust(top=0.931, bottom=0.084, left=0.101,
                            right=0.955, hspace=0.141, wspace=0.2)
        a1, a2  = fig.subplots(2,1)

        # Plot models:
        a1.plot(mag['time'], mag['dBh'],   lw=2)
        a2.plot(mag['time'], mag['dBdth'], lw=2)

        # Plot observations:
        if do_obs:
            if s in data:
                a1.plot(data['time'], data[s+'_H' ], 'k', alpha=.7)
                a2.plot(data['time'], data[s+'_dH'], 'k', alpha=.7)
        
        # Customize axes:
        applySmartTimeTicks(a1, mag['time'], dolabel=False)
        applySmartTimeTicks(a2, mag['time'], dolabel=True)
        a1.set_ylabel(r'$\Delta B_H$ ($nT$)', size=20)
        a2.set_ylabel(r'$|dB/dt|_H$ ($nT/s$)', size=20)

        # Add title:
        a1.set_title('Station {0} on {1:%Y-%m-%d}'.format(
            s, mag['time'][0]), size=20)
        
        # Add legend:
        if do_obs:
            a1.legend( ['Observations', 'SWMF'], loc='best')

        # Save figure:
        fig.savefig(path+'{}.png'.format(s))

        plt.close('all')

def plot_iono():
    '''
    Create a plot to summarize ionospheric activity.
    Requires total current calculation!
    '''

    from pickle import load
    
    # Some constants:
    bzcolor = '#3333CC'
    bycolor = '#ff9900'
    pdcolor = '#CC3300'
    north, south = '#0064b1', '#f58025'
    
    # Find log files:
    log_path = glob('./IE/IE*.log')
    ind_path = glob('./GM/geoind*.log')

    if not(log_path) or not(ind_path):
        raise ValueError('Cannot find log or geoindex file in IE/ and GM/')

    # Open currents:
    # time, nUp, nDown, nAll, nPhiMax, nPhiMin, sUp, sDown, sAll, sPhiMax, sPhiMin
    current_file=glob('./*.pkl')[0] if not args.currents else args.currents
    with open(current_file, 'rb') as f:
        curr =  load(f) 
    
    # Get that IMF file:
    imfpath = glob('./imf*.dat')[0] if not args.imffile else args.imffile
    imf = ImfInput(imfpath)

    # Open index/log files:
    log = bats.BatsLog(log_path[0])
    ind = bats.GeoIndexFile(ind_path[0])

    # Create figure and axes:
    f1 = plt.figure( figsize=(8.5,11) )
    a1, a2, a3, a4, a5 = f1.subplots(5,1, sharex=True)

    # CPCP:
    a1.plot(log['time'], log['cpcpn'], c=north, lw=2., label='Northern Hemi.')
    a1.plot(log['time'], log['cpcps'], c=south, lw=2., label='Southern Hemi.')
    a1.set_ylabel('CPCP ($kV$)')
    a1.legend(loc='best')

    # Total current:
    a2.plot(curr[0], curr[3], c=north, lw=2.)
    a2.plot(curr[0], curr[8], c=south, lw=2.)
    a2.set_ylabel('Total J$_{\parallel}$ ($MA$)')
    
    # Peak azimuthal:
    a3.plot(curr[0], curr[4]/1000., c=north, lw=2.)
    a3.plot(curr[0], curr[5]/1000., c=north, lw=2.)
    a3.plot(curr[0], curr[9]/1000., c=south, lw=2.)
    a3.plot(curr[0], curr[10]/1000.,c=south, lw=2.)
    a3.set_ylabel('Max/Min J$_{\\Phi}$ ($mA/m$)')

    # Add IMF data:
    a4.plot(imf['time'], imf['bz'], color=bzcolor, lw=2)
    a4.plot(imf['time'], imf['by'], color=bycolor, lw=2)
    a4.hlines(0, ind['time'][0], ind['time'][-1], color='k', linestyles='dashed')
    a4.legend(['B$_Z$', 'B$_Y$'], loc='best')
    a4.set_ylabel('IMF ($nT$)')

    # Add Pdyn:
    a5.plot(imf['time'], imf['pram'], color=pdcolor)
    a5.grid(False, axis='y')
    a5.set_ylabel('P$_{dyn}$ ($nPa$)')

    # Set axes ticks and x-labels:
    applySmartTimeTicks(a5, ind['time'], dolabel=True)
    for a in [a1, a2, a3, a4]:
        a.set_xlabel('')
        #a.set_xticklabels('')

    a1.set_title('Ionospheric Summary: {:%Y-%m-%d}'.format(
        log['time'][0]), size=20)
    plt.tight_layout()

    f1.savefig('./iono_summary.png')
    
def plot_indexes():
    # Some constants:
    bzcolor = '#3333CC'
    bycolor = '#ff9900'
    pdcolor = '#CC3300'
    
    # Find log files:
    log_path = glob('./GM/log*.log')
    ind_path = glob('./GM/geoind*.log')

    # Get that IMF file:
    imfpath = glob('./imf*.dat')[0] if not args.imffile else args.imffile
    imf = ImfInput(imfpath)

    if not(log_path) or not(ind_path):
        raise ValueError('Cannot find log or geoindex file in GM/')
    
    # Open index/log files:
    log = bats.BatsLog(log_path[0])
    ind = bats.GeoIndexFile(ind_path[0])

    # Create figure and axes:
    f1 = plt.figure( figsize=(8.5,11) )
    a1, a2, a3, a4, a5 = f1.subplots(5,1, sharex=True)
    
    # Set line style for observations:
    obs_kwargs={'c':'k', 'ls':'--', 'alpha':.7}
    
    # Dst:
    plot1=log.add_dst_quicklook(target=a1,plot_sym=args.obs,
                                obs_kwargs=obs_kwargs)
    
    # Kp:  
    plot2=ind.add_kp_quicklook(target=a2, loc=412, plot_obs=args.obs)
    
    # Ae:
    plot3=ind.add_ae_quicklook(target=a3, plot_obs=args.obs, val='AU',
                               obs_kwargs=obs_kwargs)
    plot3=ind.add_ae_quicklook(target=a3, plot_obs=args.obs, val='AL',
                               obs_kwargs=obs_kwargs,
                               c=a3.get_lines()[0].get_color())
    # Replace legend with something better.
    a3.set_ylabel('AU/AL ($nT$)')
    a3.legend( a3.lines[:2], ['Model AU/AL', 'Obs. AU/AL'], loc='best')
    
    # Add IMF data:
    a4.plot(imf['time'], imf['bz'], color=bzcolor, lw=1.25)
    a4.plot(imf['time'], imf['by'], color=bycolor, lw=1.25)
    a4.hlines(0, ind['time'][0], ind['time'][-1], color='k', linestyles='dashed')
    a4.legend(['B$_Z$', 'B$_Y$'], loc='best')
    a4.set_ylabel('IMF ($nT$)')

    # Add Pdyn:
    a5.plot(imf['time'], imf['pram'], color=pdcolor)
    a5.grid(False, axis='y')
    a5.set_ylabel('P$_{dyn}$ ($nPa$)')

    # Set axes ticks and x-labels:
    applySmartTimeTicks(a5, ind['time'], dolabel=True)
    for a in [a1, a2, a3, a4]:
        a.set_xlabel('')
        #a.set_xticklabels('')

    a1.set_title('Geomagnetic Indices: {:%Y-%m-%d}'.format(
        log['time'][0]), size=20)
    plt.tight_layout()

    f1.savefig('./geoindexes.png')


########################### MAIN SCRIPT #########################
if __name__ == '__main__':

    print('Plotting indexes...')
    plot_indexes()
    print('Plotting iono summary...')
    plot_iono()
    if(args.mags):
        print('Plotting magnetometers...')
        print('\tObservation file = {}'.format(args.supermag))
        plot_mags()
