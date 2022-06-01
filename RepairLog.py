#!/usr/bin/env python3
'''
Repair logfiles with multiple entries with the same timestamp.

When the SWMF has a non-standard restart (e.g., due to a crash or interruption
in the computing environment), any log file-type output will be fractured into
multiple parts. If the SWMF is restarted at a point in the simulation that
was already completed, the log file parts may have redundant or overlapping
entries. While concatenation scripts that properly handle this situation
exist (e.g., CatLog.py in this repository), some do not- merely appending
the results and keeping the overlapping time periods.

This script scans an SWMF log file (including magnetometer files, virtual sats,
and others) for overlapping time periods and removes them. The original file
will be overwritten unless the debug flag is active.

'''

import shutil
from argparse import ArgumentParser, RawDescriptionHelpFormatter


def get_time(line, index, debug=False):
    '''
    From a string entry, return the "time" the file was written.  This is
    done by splitting the line and taking all items corresponding to the
    indexes within the input list "index", concatenating them, and
    converting the resulting line into a float.  Many preceeding zeros are
    added to each entry to compensate for frequenty chances in within
    log files.
    '''

    parts = line.split()
    keep = ''

    for i in index:
        keep += '{:0>2}'.format(parts[i])

    if debug:
        print('TIME CHECK DEBUG:')
        print('Input Line="{}"'.format(line))
        print('Reducing to {}'.format(keep))

    return int(keep)


# Create argument parser & set up arguments:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-o", "--outfile", default=None,
                    help="Rather than append to first file, create new file " +
                    "to store output.")
parser.add_argument("files", nargs='+', help="Files to convert.  Can be " +
                    "explicit files or a unix wildcard.")
parser.add_argument("--debug", action="store_true",
                    help="Print debug information, keep original file.")
parser.add_argument("-k", "--keep", action="store_true",
                    help="Keep original version of files with '_orig' " +
                    "appended to the file name.")

# Handle arguments, noting that argparse expands linux wildcards.
args = parser.parse_args()

for filename in args.files:
    print(f"Repairing {filename}...")

    # If "keep"-ing originals, copy them now.
    shutil.copy(filename, filename + '_orig')

    # Slurp entire file contents.
    with open(filename, 'r') as f:
        raw = f.readlines()

    # Load and store header:
    info = raw.pop(0)  # garbage
    head = raw.pop(0)  # Header.
    nbytes = len(raw[0])

    # Using header info, create list of indexes corresponding to time.
    time_locs = []  # List of indices
    time_vars = []  # List of variable names (for debugging).

    # Desired time variable names in order:
    search_names = ['year', 'yyyy', 'yy', 'yr', 'doy', 'mo', 'month', 'mm',
                    'day', 'dy', 'hour', 'hr', 'hh', 'mm', 'mn', 'min', 'ss',
                    'sec', 'sc']
    iter_names = ['iter', 'it', 'nstep']

    # Search for time tags:
    for i, part in enumerate(head.split()):
        for s in search_names:
            if s == part.lower():
                time_locs.append(i)
                time_vars.append(s)
                break

    # If no time tags are found, try iterations:
    if not time_locs:
        for i, part in enumerate(head.split()):
            for s in iter_names:
                if s == part.lower():
                    time_locs.append(i)
                    time_vars.append(s)
                    break
            if time_locs:
                break  # Only want a single iteration tag.

    # If nothing was found still, default to first column:
    if not time_locs:
        time_locs.append(0)
        time_vars.append('Default (none found)')

    if args.debug:
        print("DEBUG:\tOpened file {}" .format(filename))
        print("\tEach line is {} characters long." .format(nbytes))
        print("\tHeader has {} entries." .format(len(head.split())))
        print("\tUsing the following columns in order for time calculation:")
        for i, s in zip(time_locs, time_vars):
            print("\t[{:02d}] {}".format(i, s))

    # Open replacement file, begin writing contents:
    last_time = 0
    with open(filename+args.debug*'_repaired', 'w') as f:
        # Write header file to new file:
        f.write(info)
        f.write(head)

        # Write data to new file. If time tag *decreases*, we have
        # overlapping entries. Skip entries until time increases again.
        for line in raw:
            time = get_time(line, time_locs, args.debug)
            if last_time > time:
                continue
            f.write(line)
            last_time = time

