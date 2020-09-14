#!/usr/bin/env python

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
parser.add_argument("-m", "--mags", default='',
                    help="Plot magnetometer observations over model results "+
                    "using values in SuperMag output file given by argument "+
                    "value (e.g, --mags 'mags.something'.")

#parser.add_argument("-obs", "--obsdst", help="Fetch and plot observed dst on "+
#                    " top of modeled values", action="store_true")
#parser.add_argument("--debug", help="Only the first plot is created.  "+
#                    "It is shown to the user and not saved.  Additional "+
#                    "information is printed to screen.", action="store_true")

# Handle arguments:
args = parser.parse_args()

from glob import glob

import numpy as np
import matplotlib.pyplot as plt

from spacepy.plot import style, applySmartTimeTicks
from spacepy.pybats import bats, ImfInput

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
    do_obs = bool(args.mags)
    if do_obs:
        data = read_supermag(args.mags)

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
        
def plot_indexes():
    # Find log files:
    log_path = glob('./GM/log*.log')
    ind_path = glob('./GM/geoind*.log')

    if not(log_path) or not(ind_path):
        raise ValueError('Cannot find log or geoindex file in GM/')
    
    # Open index/log files:
    log = bats.BatsLog(log_path[0])
    ind = bats.GeoIndexFile(ind_path[0])
    
    f1 = plt.figure( figsize=(10,10) )

    # Set line style for observations:
    obs_kwargs={'c':'k', 'ls':'--', 'alpha':.7}
    
    # Dst:
    plot1=log.add_dst_quicklook(target=f1, loc=311,
                                plot_sym=True, obs_kwargs=obs_kwargs)
    
    # Kp:  TEMPORARILY DISABLED.
    #plot2=ind.add_kp_quicklook( target=f1, loc=412, plot_obs=True)
    
    # Ae:
    plot3=ind.add_ae_quicklook( target=f1, loc=312, plot_obs=True, val='AE',
                                obs_kwargs=obs_kwargs)
    plot4=ind.add_ae_quicklook( target=f1, loc=313, plot_obs=True, val='AL',
                                obs_kwargs=obs_kwargs)

    for p in [plot1, plot3]:
        p[1].set_xlabel('')
        p[1].set_xticklabels('')

    plot1[1].set_title('Geomagnetic Indices: {:%Y-%m-%d}'.format(
        log['time'][0]), size=20)
    plt.tight_layout()

    f1.savefig('./geoindexes.png')


########################### MAIN SCRIPT #########################
if __name__ == '__main__':

    print('Plotting indexes...')
    plot_indexes()
    print('Plotting magnetometers...')
    if(args.mags):
        print('\tObservation file = {}'.format(args.mags))
    plot_mags()
