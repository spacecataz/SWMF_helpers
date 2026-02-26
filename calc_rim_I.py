#!/usr/bin/env python3

'''
From a series of Ridley_serial 2D output files, calculate the integrated
radial current density and max/min azimuthal current density for
both hemispheres and save the resulting timeseries data into a new file.

The script will look for files named "it??????_??????_???.idl" in three
separate places in the following order:
   1) ./
   2) ./IE/
   3) ./IE/ionosphere/
Once files are found, only those will be used in the calculation.  This means
that you can run this script from the top-level of an SWMF run directory.

The output file can be of two formats: simple ASCII or a Python Pickle.
The pickle will contain a dictionary of numpy arrays: an array of datetimes,
and then numpy arrays of total FACs, upward FACs, downward FACs,
maximum azimuthal current, minimum azimuthal current, total azimuthal
current, dayside FACs, nightside FACs, dayside azimuthal current, and
nightside azimuthal current, all in the northern
hemisphere; these values are then repeated for the southern
hemisphere.  The ASCII files have a self-descriptive header.

The variables are stored in the pickle-formatted files with these keys:
time, nUp, nDown, nAll, nPhiMax, nPhiMin, sUp, sDown, sAll, sPhiMax, sPhiMin,
n_J, s_J, ndayI, nnightI, sdayI, snightI, ndayJ, nnightJ, sdayJ, snightJ
...where "n" or "s" represents northern or southern hemisphere and "up",
"down", and "all" is the FAC direction ("all" is net). "day" and "night"
indicate current on either the positive or negative X side of the hemisphere,
respectively. "I" is integrated field-aligned current ("jr") and "J" is
integrated azimuthal current ("jphi").

Output units are in Mega-Amps (MA). (...at least for iFACs. Need to verify for
azimuthal current.)

Variable Summary
----------------

Each has a 'n' or 's' prefix, specifying the hemisphere.

| Name           | Value                                                      |
|----------------|------------------------------------------------------------|
| time           | Vector of datetimes                                        |
| Up, Down, All  | Integrated FAC; 'All' is net current (up + down ~= 0)      |
| PhiMax, PhiMin | Max (eastward) and min (westward) electrojet (J_phi) value |
| dayI, nightI   | Integrated absolute value of FAC per region                |
| dayJ, nightJ   | Integrated absolute value of J_phi per region              |

'''

# Load relevant modules
from glob import glob
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as np

from spacepy.pybats.rim import Iono

# Start by setting up and parsing arguments.
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('--outfile', '-o', default='./current.txt',
                    help='Set the name of the output file.  Allowable ' +
                    'suffixes are  "*.txt" or "*.pkl" for ascii or Python ' +
                    'pickle, respectively.  Default is "./current.txt"')
parser.add_argument('--debug', '-d', default=False, action='store_true',
                    help='Turn on verbose debug information.')

# Handle arguments:
args = parser.parse_args()

# Check output file name:
extension = args.outfile.split('.')[-1]
if extension != 'pkl' and extension != 'txt':
    raise ValueError('Unrecognized file format, use .pkl or .txt only.')

# Get list of files:
file_list = glob('./it??????_??????_???.idl*')
if not file_list:
    file_list = glob('./IE/it??????_??????_???.idl*')
if not file_list:
    file_list = glob('./IE/ionosphere/it??????_??????_???.idl*')

# Sort the files (thanks, OSX):
file_list.sort()

if args.debug:
    print(f"Found {len(file_list)} total files:")
    for f in file_list:
        print(f"\t{f}")

# Prepare numpy arrays to hold results:
nFiles = len(file_list)
time = np.zeros(nFiles, dtype=object)
nUp, nDown, nAll = np.zeros(nFiles), np.zeros(nFiles), np.zeros(nFiles)
sUp, sDown, sAll = np.zeros(nFiles), np.zeros(nFiles), np.zeros(nFiles)

nPhiMax, nPhiMin = np.zeros(nFiles), np.zeros(nFiles)
sPhiMax, sPhiMin = np.zeros(nFiles), np.zeros(nFiles)

n_J, s_J = np.zeros(nFiles), np.zeros(nFiles)

ndayI, nnightI = np.zeros(nFiles), np.zeros(nFiles)
sdayI, snightI = np.zeros(nFiles), np.zeros(nFiles)

ndayJ, nnightJ = np.zeros(nFiles), np.zeros(nFiles)
sdayJ, snightJ = np.zeros(nFiles), np.zeros(nFiles)


# Create function that integrates horizontal current, day/night FACs:
def integrate_more(ie):
    '''
    Takes a spacepy.pybats.rim Iono object that has already
    used the method ie.calc_j()

    Returns the value of total integrated horizontal current,
    the integrated dayside FACs, the integrated nightside FACs,
    the integrated dayside jphi, and the integrated nightside
    jphi.

    Modeled strongly after the integration in ie.calc_I()
    '''

    # Calculate some physically meaningful values/units
    units = 1E-6*1E-6  # micro amps to amps, amps to MegaAmps
    R = (6371.0+110.0)*1000.0  # Radius of Earth + iono altitude
    dTheta = np.pi*ie.dlat/180.
    dPhi = np.pi*ie.dlon/180.

    def make_integrand(values, colat):
        '''For integrating a value in IE over the whole hemisphere.'''
        return values*np.sin(colat)*dTheta*dPhi

    ieshape = ie['n_phi'].shape

    for hemi in 'ns':
        colat = ie[hemi+'_theta']*np.pi/180.  # for integration later

        has_jphi = 'n_jx' in ie
        has_pos = 'n_x' in ie

        currents = ['jr'] + ['jphi'] * has_jphi

        # Split into dayside and nightside:
        day = ie[hemi+'_x'] > 0. if has_pos else np.zeros(ieshape, dtype=bool)
        night = ie[hemi+'_x'] < 0 if has_pos else np.zeros(ieshape, dtype=bool)

        for current in currents:
            day_j = np.copy(ie[hemi+'_'+current])
            night_j = np.copy(ie[hemi+'_'+current])
            day_j[night] = 0.
            night_j[day] = 0.

            # Get locations of "up" and "down" (really, positive and negative):
            loc_up = ie[hemi+'_'+current] > 0
            loc_do = ie[hemi+'_'+current] < 0

            # Integrate day/night currents ONLY if position available.
            for val, name in zip([day_j, night_j], ['day', 'night']):
                integrand = make_integrand(val, colat)
                key = hemi+'_'+name+'_'+current  # for convenience

                # Integrate "up" and "down" day or night current:
                ie[key+'up'] = units*R**2 * np.sum(integrand[loc_up])
                ie[key+'down'] = units*R**2 * np.sum(integrand[loc_do])

                # Get total day or night current:
                ie[key] = (np.abs(ie[key+'up']) + np.abs(ie[key+'down']))/2

        # Integrate total jphi
        if has_jphi:
            loc_up = ie[hemi+'_jphi'] > 0
            loc_do = ie[hemi+'_jphi'] < 0

            integrand = make_integrand(ie[hemi+'_jphi'], colat)
            ie[hemi+'_up_jphi'] = units*R**2 * np.sum(integrand[loc_up])
            ie[hemi+'_down_jphi'] = units*R**2 * np.sum(integrand[loc_do])

            ie[hemi+'_J'] = (ie[hemi+'_up_jphi']
                             + np.abs(ie[hemi+'_down_jphi']))/2
        else:
            ie[hemi+'_jphi'] = np.zeros(ieshape)
            ie[hemi+'_day_jphi'], ie[hemi+'_night_jphi'] = 0, 0
            ie[hemi+'_J'] = 0.0

    return ie


# Load up variables by looping over all files:
for i, f in enumerate(file_list):
    ie = Iono(f)  # open file.

    ie.calc_I()   # integrate current.
    if 'n_jx' in ie:
        ie.calc_j()   # get azimuthal current.

    ie = integrate_more(ie)  # integrate horizontal current and day/night FACs

    # Store time and north/south hemi integrated currents:
    time[i] = ie.attrs['time']
    nUp[i], nDown[i], nAll[i] = ie['n_Iup'], ie['n_Idown'], ie['n_I']
    sUp[i], sDown[i], sAll[i] = ie['s_Iup'], ie['s_Idown'], ie['s_I']

    # Store maximum/minimum azimuthal currents:
    nPhiMax[i], nPhiMin[i] = ie['n_jphi'].max(), ie['n_jphi'].min()
    sPhiMax[i], sPhiMin[i] = ie['s_jphi'].max(), ie['s_jphi'].min()

    # Store horizontal currents:
    n_J[i], s_J[i] = ie['n_J'], ie['s_J']

    # Store day/night FACs:
    ndayI[i], nnightI[i] = ie['n_day_jr'], ie['n_night_jr']
    sdayI[i], snightI[i] = ie['s_day_jr'], ie['s_night_jr']

    # Store day/night Jphi:
    ndayJ[i], nnightJ[i] = ie['n_day_jphi'], ie['n_night_jphi']
    sdayJ[i], snightJ[i] = ie['s_day_jphi'], ie['s_night_jphi']


# Now, save as a file:
if args.outfile.split('.')[-1] == 'pkl':
    import pickle
    f = open(args.outfile, 'wb')
    data = {'time': time,
            'nUp': nUp, 'nDown': nDown, 'nAll': nAll,
            'sUp': sUp, 'sDown': sDown, 'sAll': sAll,
            'nPhiMax': nPhiMax, 'nPhiMin': nPhiMin,
            'sPhiMax': sPhiMax, 'sPhiMin': sPhiMin,
            'n_J': n_J, 's_J': s_J, 'ndayI': ndayI, 'nnightI': nnightI,
            'sdayI': sdayI, 'snightI': snightI, 'ndayJ': ndayJ,
            'nnightJ': nnightJ, 'sdayJ': sdayJ, 'snightJ': snightJ}
    pickle.dump(data, f)
    f.close()
else:
    f = open(args.outfile, 'w')
    # Header:
    f.write('Time\tI_Up_North(MA)\tI_Down_North(MA)\tI_Total_North(MA)')
    f.write('\tJphi_Max_North(A/m)\tJphi_Min_North(A/m)')
    f.write('\tI_Up_South(MA)\tI_Down_South(MA)\tI_Total_South(MA)')
    f.write('\tJphi_Max_South(A/m)\tJphi_Min_South(A/m)')
    f.write('\tTotal_J_North(A/m)\tTotal_J_South(A/m)')
    f.write('\tI_Day_North(A/m)\tI_Night_North(A/m)')
    f.write('\tI_Day_South(A/m)\tI_Night_South(A/m)')
    f.write('\tJ_Day_North(A/m)\tJ_Night_North(A/m)')
    f.write('\tJ_Day_South(A/m)\tJ_Night_South(A/m)\n')

    for i, t in enumerate(time):
        tnow = t.isoformat()
        f.write('{}\t'.format(tnow))
        f.write('{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.5f}\t'.format(
            nUp[i], nDown[i], nAll[i], nPhiMax[i], nPhiMin[i]))
        f.write('{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.5f}\t'.format(
            sUp[i], sDown[i], sAll[i], sPhiMax[i], sPhiMin[i]))
        f.write('{:12.5f}\t{:12.5f}\t'.format(n_J[i], s_J[i]))
        f.write('{:12.5f}\t{:12.5f}\t'.format(ndayI[i], nnightI[i]))
        f.write('{:12.5f}\t{:12.5f}\t'.format(sdayI[i], snightI[i]))
        f.write('{:12.5f}\t{:12.5f}\t'.format(ndayJ[i], nnightJ[i]))
        f.write('{:12.5f}\t{:12.5f}\n'.format(sdayJ[i], snightJ[i]))
    f.close()
