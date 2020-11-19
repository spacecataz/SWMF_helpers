#!/usr/bin/env python

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
The pickle will contain a dictionary of 7 numpy arrays: an array of datetimes,
and then numpy arrays of total current, upward current, downward current,
maximum azimuthal current, and minimum azimuthal current, all in the northern 
hemisphere; these values are then repeated for the southern
hemisphere.  The ASCII files have a self-descriptive header.

The variables are stored in the pickle-formatted files with these keys:
time, nUp, nDown, nAll, nPhiMax, nPhiMin, sUp, sDown, sAll, sPhiMax, sPhiMin
...where "n" or "s" represents northern or southern hemisphere and "up", 
"down", and "all" is the FAC direction ("all" is net.)

Output units are in Mega-Amps (MA).

Failure may occur due to formatting errors in the Ridley_serial output files.
This can be addressed effectively with the "--fix" argument (see below).
'''

# Start by setting up and parsing arguments.
from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('--outfile', '-o', default='./current.txt', help=
                    'Set the name of the output file.  Allowable suffixes '+
                    'are  "*.txt" or "*.pkl" for ascii or Python pickle, ' +
                    'respectively.  Default is "./current.txt"')
parser.add_argument('--fix', '-f', default=False, action='store_true', help=
                    'Ridley_serial 2D output has poor format specification.'+
                    '  This can be fixed by SpacePy.  If this program fails '+
                    'while reading a 2D file, use this option to repair all '+
                    'broken files [warning: slow].')
parser.add_argument('--debug', '-d', default=False, action='store_true', help=
                    'Turn on verbose debug information.')

# Handle arguments:
args = parser.parse_args()

# Load relevant modules
from glob import glob
import numpy as np
from spacepy.pybats.rim import Iono, fix_format

# Check output file name:
if args.outfile.split('.')[-1]!='pkl' and args.outfile.split('.')[-1]!='txt':
    raise ValueError('Unrecognized file format, use .pkl or .txt only.')

# Get list of files:
file_list = glob('./it??????_??????_???.idl')
if not file_list:
    file_list=glob('./IE/it??????_??????_???.idl')
if not file_list:
    file_list=glob('./IE/ionosphere/it??????_??????_???.idl')

# Sort the files (thanks, OSX):
file_list.sort()
    
# Fix files if requested:
if args.fix: 
    for f in file_list:
        fix_format(f)

# Prepare numpy arrays to hold results:
nFiles = len(file_list)
time = np.zeros(nFiles, dtype=object)
nUp, nDown, nAll = np.zeros(nFiles), np.zeros(nFiles), np.zeros(nFiles)
sUp, sDown, sAll = np.zeros(nFiles), np.zeros(nFiles), np.zeros(nFiles)

nPhiMax, nPhiMin = np.zeros(nFiles), np.zeros(nFiles)
sPhiMax, sPhiMin = np.zeros(nFiles), np.zeros(nFiles)

# Load up variables by looping over all files:
for i, f in enumerate(file_list):
    ie = Iono(f)  #open file.
    ie.calc_I()   #integrate current.
    ie.calc_j()   #get azimuthal current.
    
    # Store time and north/south hemi integrated currents:
    time[i] = ie.attrs['time']
    nUp[i], nDown[i], nAll[i] = ie['n_Iup'], ie['n_Idown'], ie['n_I']
    sUp[i], sDown[i], sAll[i] = ie['s_Iup'], ie['s_Idown'], ie['s_I']

    # Store maximum/minimum azimuthal currents:
    nPhiMax[i], nPhiMin[i] = ie['n_jphi'].max(), ie['n_jphi'].min()
    sPhiMax[i], sPhiMin[i] = ie['s_jphi'].max(), ie['s_jphi'].min()
    
# Now, save as a file:
if args.outfile.split('.')[-1] == 'pkl':
    import pickle
    f = open(args.outfile, 'wb')
    data = {'time':time,
            'nUp':nUp,'nDown':nDown,'nAll':nAll,'nPhiMax':nPhiMax,'nPhiMin':nPhiMin,
            'sUp':sUp,'sDown':sDown,'sAll':sAll,'sPhiMax':sPhiMax,'sPhiMin':sPhiMin}
    pickle.dump(data, f)
    f.close()
else:
    f = open(args.outfile, 'w')
    # Header:
    f.write('Time\tI_Up_North(MA)\tI_Down_North(MA)\tI_Total_North(MA)')
    f.write(    '\tJphi_Max_North(A/m)\tJphi_Min_North(A/m)')
    f.write(    '\tI_Up_South(MA)\tI_Down_South(MA)\tI_Total_South(MA)')
    f.write(    '\tJphi_Max_South(A/m)\tJphi_Min_South(A/m)')

    for i, t in enumerate(time):
        tnow = t.isoformat()
        f.write('{}\t{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.5f}\t'.format(
            tnow, nUp[i], nDown[i], nAll[i], nPhiMax[i], nPhiMin[i]))
        f.write('{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.5f}\n'.format(
            sUp[i], sDown[i], sAll[i], sPhiMax[i], sPhiMin[i]))
    f.close()
