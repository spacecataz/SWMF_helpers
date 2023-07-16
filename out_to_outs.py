#!/usr/bin/env python3

'''
This script will identify sets of files to be converted from individual
".out" files to a single ".outs" file.

The SWMF can produce ".out" files, which represent output from a single point
in time. These can be concatenated into ".outs" files (note the "s") which
contain all data from all individual ".out" files. While it is convenient to
produce .outs files from within an existing run directory, it is less easy to
do so after the data is moved to a different location.

This script quickly identifies different groups of output files within an
SWMF output directory and concatenates output into ".outs".
'''

import os
import re
import sys
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument('--remove', '-rm', help="Remove .out files after " +
                    'concatenation leaving only the .outs.',
                    default=True, action='store_true')
parser.add_argument('--debug', '-d', default=False, action='store_true',
                    help='Turn on verbose debug information.')
# Handle arguments:
args = parser.parse_args()

# Start by collecting file types that can be concatenated.
out_raw = glob('GM/*.out') + glob('GM/IO2/*.out') + glob('*.out')

# Loop over files and determine what groups exist.
ftypes = {}
print('Searching for file groups: ')
for f in out_raw:
    x = re.search('(.+\_[ent])\d+.*.out', f)
    if x:
        if x[1] not in ftypes:
            ftypes[x[1]] = [f]
        else:
            ftypes[x[1]] += [f]

# Sort file lists; report back to user.
for f in ftypes:
    ftypes[f].sort()
    print(f"\t{f}\t({len(ftypes[f])} files found)")

# Ask to continue before proceeding
if args.debug:
    response = input('Continue [Y/n]?')
    if response == 'n':
        sys.exit()


# Now, concatenate.
for ftype in ftypes:
    print(f"Working on group {ftype}")

    # Get first and last time stamp values.
    nprefix = len(ftype)
    stamp1 = ftypes[ftype][0][nprefix:-4]
    stamp2 = ftypes[ftype][-1][nprefix:-4]
    if args.debug:
        print('First and last timestamp:')
        print(f'\t{stamp1}\n\t{stamp2}')

    # Create new file; concatenate each file into it.
    outfile = f"{ftype}{stamp1}_{stamp2}.outs"
    print(f"Creating new file {outfile}")
    with open(outfile, 'wb') as out:
        # Loop over all files in type
        for filenow in ftypes[ftype]:
            # Dump into output file
            if args.debug:
                print(f"\tConcatenating {filenow}...")
            with open(filenow, 'rb') as f:
                out.write(f.read())
        # Only remove files once writing of outs is complete:
        if args.remove:
            print('\tRemoving individual files...')
            for filenow in ftypes[ftype]:
                os.remove(filenow)
