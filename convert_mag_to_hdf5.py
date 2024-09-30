#!/usr/bin/env python

'''
Convert a magnetomter `.outs` file group (note the "s" in "outs"!) into a
simplified HDF file for use in the Solar Tsunamis project.
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import h5py
from matplotlib.dates import date2num, get_epoch
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
args = parser.parse_args()

if args.verbose:
    print(f'Reading mag grid file: {args.file}')
if 'mags' not in globals():
    mags = MagGridFile(args.file)

# Collect info about file:
ntime, nlat, nlon = mags.attrs['nframe'], mags['Lat'].size, mags['Lon'].size

# Create HDF5 file:
out = h5py.File(args.outfile + '.h5', 'w')

# Create datasets for position:
for v, size in zip(['lon', 'lat'], [nlon, nlat]):
    set = out.create_dataset(v, [size])
    set[:] = mags[v.capitalize()]

# Create datasets for perturbation and populate:
for v in ['dBn', 'dBe', 'dBd']:
    out.create_dataset(v, [nlon, nlat, ntime])

# Create time array:
set = out.create_dataset('time', [ntime])
set[:] = date2num(mags.attrs['times'])
set = out.create_dataset('ref_epoch', shape=1, dtype=h5py.string_dtype())
set = get_epoch()

# Populate and close:
for i in range(ntime):
    mags.switch_frame(i)
    for v in ['dBn', 'dBe', 'dBd']:
        out[v][:, :, i] = mags[v]

out.close()