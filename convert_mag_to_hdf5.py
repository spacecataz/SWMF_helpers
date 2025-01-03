#!/usr/bin/env python3

'''
Convert a magnetomter `.outs` file group (note the "s" in "outs"!) into a
simplified set of HDF files for use in the Solar Tsunamis project.

Note that multiple HDF files will be created to break up the data into
manageable chunks. The `chunk` argument controls this behavior.

Files will be saved as {prefix}_n{filenumber}.hdf5
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import h5py
from numpy import ceil
from spacepy.pybats.bats import MagGridFile

# Start by configuring the argparser:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("file", type=str,
                    help='A ".outs" file to convert.')
parser.add_argument("-v", "--verbose", default=False, action='store_true',
                    help="Turn on verbose output mode.")
parser.add_argument("--debug", default=False, action='store_true',
                    help="Turn on debugging mode: create a small file for " +
                    "testing & validation purposes.")
parser.add_argument("-p", "--prefix", default='maggrid', help="Set output " +
                    "file name without file extension. Defaults to 'maggrid'")
parser.add_argument("-c", "--chunk", default=60, help="Set the size of file" +
                    " time chunks, in minutes. 60 means 1 file will contain " +
                    "60 minutes of data. If negative, results are not " +
                    "chunked.", type=int)
parser.add_argument("--nvar", "-n",
                    choices=['min', 'max', 'med'], default='min',
                    help="Set number of variables to save. 'min' will save " +
                    "only 1 component (dBd); 'med' saves all 3 components, " +
                    "'max' saves all componets broken by contribution.")
args = parser.parse_args()

# turn on Verbose in debug mode:
if args.debug:
    args.verbose = True

if args.verbose:
    print(f'Reading mag grid file: {args.file}')
if 'mags' not in globals():
    mags = MagGridFile(args.file)

# Collect info about file:
ntime, nlat, nlon = mags.attrs['nframe'], mags['Lat'].size, mags['Lon'].size
if args.verbose:
    print(f"nTimes={ntime}, nLat x nLon = {nlat} x {nlon}")

# Create short file for debug mode:
if args.debug:
    print('DEBUG: Shortening file to 10 time entries.')
    ntime = 10

# Determine the number of files required to write based on chunk size:
dtime = (mags.attrs['times'][1] - mags.attrs['times'][0]).total_seconds()

# Get number of entries per chunk. Only 1 file if args.chunk is negative.
strsize = 'ALL' if args.chunk < 0 else args.chunk
chunksize = ntime if args.chunk < 0 else args.chunk / (dtime / 60.)

# Get number of files: total time in hours / chunk size.
nfiles = int(ceil(ntime / chunksize))

if args.verbose:
    print(f"Each entry in the input file covers {dtime}s.")
    print(f"Each output file will contain {strsize} hours.")
    print(f"This equates to {chunksize} entries per HDF and " +
          f"{nfiles} total HDF files.")

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

if args.verbose:
    print('Saving the following variables to file:')
    for v in savevars:
        print(f'\t{v}')

# Generate time as floating point:
starttime = mags.attrs['times'][0]
time = [(t - starttime).total_seconds() for t in mags.attrs['times']]

# Loop through files, grabbing values we want.
for ifile in range(nfiles):
    # Create HDF5 file:
    if args.verbose:
        print(f'Creating HDF5 file {args.prefix}_n{ifile:05d}.h5')

    # Get start/end time for values to store in this file.
    istart = int(ifile * chunksize)
    iend = min(ntime, int((ifile + 1) * chunksize))

    # Last file may not be a complete chunk. Address this.
    chunknow = iend - istart

    # Check'em!
    if args.verbose:
        print(f'Working on entries [{istart}:{iend}) (chunksize={chunknow})')

    with h5py.File(f"{args.prefix}_n{ifile:05d}.h5", 'w') as out:
        # Create datasets for position:
        for v, size in zip(['lon', 'lat'], [nlon, nlat]):
            set = out.create_dataset(v, [size])
            set[:] = mags[v.capitalize()]

        # Create datasets for perturbation and populate: [ntime, nlon, nlat]
        for v in savevars:
            out.create_dataset(v, [chunknow, nlon, nlat])

        # Create time array:
        set = out.create_dataset('time', [chunknow])
        set[:] = time[istart:iend]
        set = out.create_dataset('ref_epoch', shape=1,
                                 dtype=h5py.string_dtype())
        set[:] = starttime.isoformat()

        # Populate and close:
        if args.verbose:
            print(f"Writing file number {ifile+1} of {nfiles}")
        for i in range(istart, iend):
            if args.verbose:
                print(f'\t{i/ntime:7.2%} complete ' +
                      f'({i+1} of {ntime} written...)')
            mags.switch_frame(i)
            for v in savevars:
                out[v][i-istart, :, :] = mags[v]

# In debug mode, verify file:
if args.debug:
    print('Verifying file contents...')
    with h5py.File(f"{args.prefix}_n{0:05d}.h5", 'r') as filein:
        if filein['dBd'].shape != (chunknow, nlon, nlat):
            print('VERIFICATION FAILED. Check output.')
        if filein['ref_epoch'][0].decode() != starttime.isoformat():
            print('VERIFICATION FAILED. Check output.')

    filein.close()