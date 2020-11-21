#!/usr/bin/env python
'''
Given a series of magnetometer grid output from the SWMF, extract the results
at certain lat/lon locations and save as a single magnetometer output 
timeseries file.

The script will look for files named "mag_grid*.out" in three
separate places in the following order:
   1) <dir>/
   2) <dir>/GM/
   3) <dir>/GM/IO2/
...where <dir> defaults to "./" but can be changed via the "--dir" flag.
Once files are found, only those will be used in the calculation.  This means
that you can run this script from the top-level of an SWMF run directory.

Output file is stored in the same directory as the source files.

EXAMPLE USAGE: 
maggrid_extract.py 120 125 15 20 (Extracts from lon=120-125, lat=15-20)
'''

# Start by setting up and parsing arguments.
from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('lons', nargs=2, type=float, help='Set the latitude ' +
                    'range over which to extract.')
parser.add_argument('lats', nargs=2, type=float, help='Set the latitude ' +
                    'range over which to extract.')
parser.add_argument('--lonskip', '-xs', default=1, type=int, help=
                    'Set the cadence with which to sample. A cadence of 2 ' +
                    'will skip every other longitude entry.')
parser.add_argument('--latskip', '-ys', default=1, type=int, help=
                    'Set the cadence with which to sample. A cadence of 2 ' +
                    'will skip every other latitude entry.')
parser.add_argument('--tstart', '-t1', default=None, help=
                    'Set the start time.  Only files with timestamp at or '+
                    'after the start time will be included.  ' +
                    'String format should be YYYY-MM-DDTHH:MN:SS')
parser.add_argument('--tstop', '-t2', default=None, help=
                    'Set the stop time.  Only files with timestamp at or '+
                    'before the stop time will be included.  ' +
                    'String format should be YYYY-MM-DDTHH:MN:SS')
parser.add_argument('--dir', default='.', help='Set the name of the run '+
                    'directory in which to search for files.  Defaults to PWD.')
parser.add_argument('--outfile', '-o', default='extracted_mags.mag', help=
                    'Set the name of the output file.  '+
                    'Defaults to extracted_mags.mag')
parser.add_argument('--info', '-i', default=False, action='store_true', help=
                    'Print grid info to stdout, quit.')
parser.add_argument('--debug', '-d', default=False, action='store_true', help=
                    'Turn on verbose debug information.')

# Handle arguments:
args = parser.parse_args()

# Load relevant modules
from glob import glob
import numpy as np
from dateutil.parser import parse
from spacepy.pybats.bats import MagGridFile
from spacepy.pybats import parse_filename_time

# Get list of files, set output directory.
file_list = glob(args.dir+'/mag_grid*.out')
outdir    = args.dir+'/'
if not file_list:
    file_list=glob(args.dir+'/GM/mag_grid*.out')
    outdir = outdir+'GM/'
if not file_list:
    file_list=glob(args.dir+'/GM/IO2/mag_grid*.out')
    outdir= outdir+'GM/IO2/'
if not file_list:
    raise ValueError(f"Could not find any MagGrid files in {args.dir}")


###### CONSTANTS ##########
varlist = 'dBn dBe dBd dBnMhd dBeMhd dBdMhd dBnFac dBeFac '+\
          'dBdFac dBnHal dBeHal dBdHal dBnPed dBePed dBdPed'

###### SET UP #######
nfiles_all = len(file_list)

out = open(outdir+args.outfile, 'w')

# Get file times:
t1_all = parse_filename_time(file_list[ 0])[-1]
t2_all = parse_filename_time(file_list[-1])[-1]

# Get time from files OR flags:
if not args.tstart:
    t1 = t1_all
else:
    t1 = parse(args.tstart)

if not args.tstop:
    t2 = t2_all
else:
    t2 = parse(args.tstop)

# Reduce file list appropriately:
if t1 != t1_all:
    for i in range(nfiles_all):
        tnow = parse_filename_time(file_list[0])[-1]
        #print(t1, tnow) #debugging line
        if tnow<t1: file_list.pop(0)
        if tnow>=t1: break
nfiles = len(file_list)
if t2 != t2_all:
    for i in range(nfiles):
        tnow = parse_filename_time(file_list[-1])[-1]
        #print(t2, tnow) #debugging line
        if tnow>t2: file_list.pop(-1)
        if tnow<=t2: break
nfiles = len(file_list)

# Open first file, get information:
mag  = MagGridFile(file_list[0])
lons = mag['Lon']
lats = mag['Lat']
dLon = lons[1]-lons[0]
dLat = lats[1]-lats[0]

if args.debug or args.info:
    print(f'Mag Grid File info:')
    print(f'\tNumber of available files = {nfiles_all}')
    print(f'\tFiles have {lons.size} X {lats.size} lon-lats')
    print(f'\tLon range is {lons[0]} -- {lons[-1]} (dLon={dLon})')
    print(f'\tLat range is {lats[0]} -- {lats[-1]} (dLat={dLat})')
    print(f'\tFirst file time = {t1_all}')
    print(f'\tLast file time  = {t2_all}')
    
if args.info: exit()

# Create indices of different lats/lons:
iall = np.array(range(lons.size))
jall = np.array(range(lats.size))

# Reduce to only positions we care about:
xloc = (lons>=args.lons[0])&(lons<=args.lons[1])
yloc = (lats>=args.lats[0])&(lats<=args.lats[1])
# Extracted indexes, lat/lons and number of extractions in each:
ilon, jlat = iall[xloc][::args.lonskip], jall[yloc][::args.latskip]
elon, elat = lons[ilon], lats[jlat]
nlon, nlat = elon.size, elat.size
nMags = nlon*nlat


if args.debug:
    print(f'DEBUG INFO:')
    print(f'\tExtracting in longitude from {elon[0]} to {elon[-1]}.')
    print(f'\tExtracting in latitude  from {elat[0]} to {elat[-1]}.')
    print(f'\tLon/Lat skip = {args.lonskip}, {args.latskip}.')
    print(f'\tExtracting {nlon} X {nlat} lons/lats.')
    print(f'\t...for a total of {nMags} extractions per file.')
    print(f'\tExtraction starts at T={t1}')
    print(f'\t      ...and ends at T={t2}.')
    print(f'\tTotal number of files to be handled = {nfiles}.')


    answer = input('CONTINUE? [y/N]:')
    if answer.lower()[0] != 'y': exit()

# Create magnetometer names.  3-letter names be damned.
# Write names to file one at a time.
out.write(f'\t{nMags} magnetometers: ')
names = np.zeros( [nlon, nlat], dtype='U11')
for i, lon in enumerate(elon):
    for j, lat in enumerate(elat):
        l_now = np.abs(lat) # Get absolute value of latitude.
        prefix = 's' if lat< 0 else 'n' # Set north/south prefix
        names[i,j] = f'mag_{lon:03.0f}_{prefix}{l_now:02.0f}'
        out.write(f' {names[i,j]}')
out.write('\n')

# Finish writing header:
out.write('nstep year mo dy hr mn sc msc station X Y Z '+varlist+'\n')

# Write magnetometer file:
for f in file_list:
    # Open file:
    mag = MagGridFile(f)

    # Loop over all lons/lats:
    for i, il in enumerate(list(ilon)):
        for j, jl in enumerate(jlat):
            
            # Write time and iteration:
            out.write('{:8d} '.format(mag.attrs['iter']))
            out.write('{:%Y %m %d %H %M %S} 000 '.format(mag.attrs['time']))

            # Write mag number:
            imag = (i+1)*(j+1)
            out.write(f'{imag:4d}')
            
            # Fake XYZ:
            out.write(3*' {:12.5E}'.format(0.0))
            
            # Write the rest of the record:
            for v in varlist.split():
                out.write(' {:12.5E}'.format(mag[v][il,jl]))
            out.write('\n')
out.close()
    
