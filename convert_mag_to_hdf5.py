#!/usr/bin/env python3

'''
Convert a magnetomter `.outs` file group (note the "s" in "outs"!) into a
simplified HDF file for use in the Solar Tsunamis project.
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import h5py
from spacepy.pybats.bats import MagGridFile

# Start by configuring the argparser:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("file", type=str,
                    help='A ".outs" file to convert.')
parser.add_argument("-v", "--verbose", default=False, action='store_true',
                    help="Turn on verbose output mode.")
parser.add_argument("--debug", default=False, action='store_true',
                    help="Turn on debugging mode.")
parser.add_argument("-o", "--outfile", default='maggrid', help="Set output " +
                    "file name without file extension Defaults to 'maggrid'")
parser.add_argument("--nvar", "-n",
                    choices=['min', 'max', 'med'], default='min',
                    help="Set number of variables to save. 'min' will save " +
                    "only 1 component (dBd); 'med' saves all 3 components, " +
                    "'max' saves all componets broken by contribution.")
args = parser.parse_args()

if args.verbose:
    print(f'Reading mag grid file: {args.file}')
if 'mags' not in globals():
    mags = MagGridFile(args.file)

# Collect info about file:
ntime, nlat, nlon = mags.attrs['nframe'], mags['Lat'].size, mags['Lon'].size

# Set variables to use:
match args.nvar:
    case 'min':
        savevars = ['dBd']
    case 'med':
        savevars = ['dBn', 'dBe', 'dBd']
    case 'max':
        savevars = ['dBn', 'dBe', 'dBd', 'dBnMhd', 'dBeMhd', 'dBdMhd',
                    'dBnFac', 'dBeFac', 'dBdFac', 'dBnHal', 'dBeHal',
                    'dBdHal', 'dBnPed', 'dBePed', 'dBdPed']

if args.debug:
    print('Saving the following variables to file:')
    for v in savevars:
        print(f'\t{v}')

# Create HDF5 file:
if args.debug:
    print(f'Creating HDF5 file {args.outfile + ".h5"}')
out = h5py.File(args.outfile + '.h5', 'w')

# Create datasets for position:
for v, size in zip(['lon', 'lat'], [nlon, nlat]):
    set = out.create_dataset(v, [size])
    set[:] = mags[v.capitalize()]

# Create datasets for perturbation and populate:
for v in savevars:
    out.create_dataset(v, [nlon, nlat, ntime])

# Generate time as floating point:
starttime = mags.attrs['times'][0]
time = [(t - starttime).total_seconds() for t in mags.attrs['times']]

# Create time array:
set = out.create_dataset('time', [ntime])
set[:] = time
set = out.create_dataset('ref_epoch', shape=1, dtype=h5py.string_dtype())
set = starttime.isoformat()

# Populate and close:
for i in range(ntime):
    mags.switch_frame(i)
    for v in savevars:
        out[v][:, :, i] = mags[v]

out.close()