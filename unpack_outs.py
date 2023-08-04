#!/usr/bin/env python3

'''
This script unpacks SWMF-formatted ".outs" files into individual ".out"
components.

The SWMF can produce ".out" files, which represent output from a single point
in time. These can be concatenated into ".outs" files (note the "s") which
contain all data from all individual ".out" files. On certain occasions, it
is desireable to unpack ".outs" files into a series of ".out" files. This
script accomplishes that.

Example usage:
unpack_outs.py ./mag_grid_n00008000_00480948.outs

'''

import os
import re
import datetime as dt
from argparse import ArgumentParser
import numpy as np
from spacepy import pybats

parser = ArgumentParser(description=__doc__)
parser.add_argument("files", nargs='+', help="*.outs files to unpack. " +
                    "Can be explicit files or a unix wildcard.")
parser.add_argument('--remove', '-rm', help="Remove original .outs file " +
                    'after concatenation leaving only the .outs.',
                    default=False, action='store_true')
parser.add_argument('--verbose', '-v', default=False, action='store_true',
                    help='Turn on verbose mode.')
parser.add_argument('--debug', '-d', default=False, action='store_true',
                    help='Turn on debug information.')

args = parser.parse_args()
if args.debug:
    args.verbose = True


def gen_filename(prefix, date=None, time=None, iter=None):
    '''
    Given a prefix, iteration/date/runtime, produce a file name.
    Prefix should be a string; iter an integer, date a datetime, and
    time should be integer seconds from start of simulation.
    '''

    name = prefix

    if type(date) is dt.datetime:
        name = name + f"_e{date:%Y%M%D-%H%M%S}-000"
    elif time is not None:
        hour = int(time/3600)
        mins = int((time-3600*hour)/60)
        secs = int(time - 3600*hour - 60*mins)
        name = name + f"_t{hour:04d}{mins:02d}{secs:02d}"
    elif iter:
        name = name + f"_n{iter:08d}"

    return name + '.out'


# Loop over each file in the argument list. Open, scan, and prep to unpack.
if args.verbose:
    print(f'Gathering information for unpacking {len(args.files)} files...')
info = {}
for f in args.files:
    # Check for correct format:
    if '.outs' not in f:
        raise ValueError(f'{f} is not an ".outs" file.')

    # Open file to get critical information:
    temp = pybats.IdlFile(f)

    # Stash relevant info:
    info[f] = temp.attrs
    info[f]['offsets'] = temp._offsets

    # Get file name prefix:
    x = re.search('(.+)\_[ent]\d+.*.outs', f)
    if not x:
        raise ValueError(f'File does not appear to be valid outs file: {f}')
    info[f]['prefix'] = x[1]

    # Close the file.
    del temp

if args.verbose:
    print('Files to be unpacked:')
    for f in args.files:
        print(f"\t{f} ({info[f]['nframe']} frames)")

if args.debug:
    print('Original files will be '
          + args.remove*'REMOVED.' + (not args.remove)*'RETAINED.')
    if input('Continue [y/N]?').lower() != 'y':
        print('Aborting.')
        exit()

# Unpack each .outs file:
for f in args.files:
    if args.verbose:
        print(f'Unpacking {f}...')

    # Open the binary file:
    start_t = dt.datetime.now()
    with open(f, 'rb') as outs:
        # Create size of each frame:
        framesize = np.zeros(info[f]['nframe'], dtype=int)
        framesize[:-1] = info[f]['offsets'][1:] - info[f]['offsets'][:-1]

        # Set the last chunk size:
        outs.seek(0, 2)
        framesize[-1] = outs.tell() - info[f]['offsets'][-1]

        # Loop over each frame:
        for i in range(info[f]['nframe']):
            # Construct output file name:
            outname = gen_filename(info[f]['prefix'], info[f]['times'][i],
                                   info[f]['runtimes'][i], i)
            if args.verbose:
                print(f'\tWriting {outname}...')

            # Fast-forward to start of frame and load it:
            outs.seek(info[f]['offsets'][i])
            raw = outs.read(framesize[i])

            # Dump to file
            with open(outname, 'wb') as out:
                out.write(raw)

    # Grab timing info:
    info[f]['t'] = (dt.datetime.now() - start_t).total_seconds()
    if args.verbose:
        print(f"\tFINISHED in {info[f]['t']}s")

if args.debug:
    print('Timing info:')
    print(f"{'File Prefix':^20s}\t{'Time (s)'}")
    for f in args.files:
        print(f"{info[f]['prefix']:>20s}\t{info[f]['t']:07.1f}")

if args.remove and not args.debug:
    for f in args.files:
        os.remove(f)
