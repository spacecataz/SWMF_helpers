#!/usr/bin/env python3

'''
For each file given in argument list, count the number of files listed in
that directory and report that to the user.  Results are sorted so that the
directory with the largest number of files is listed first.

Unix wildcards are allowed.

Usage:
    countfiles.py [opts] dir1 dir2 ...

Options:
    -h or -help: print this help
'''

import os
from sys import argv
from glob import glob

dirs = []

# Parse arguments.
for arg in argv[1:]:
    if arg.lower()[:2] == '-h':
        print(__doc__)
        exit()
    else: 
        dirs+= glob(arg)

# Count files:
results = {}
total   = 0
for d in dirs:
    results[d] = sum([len(files) for x1, x2, files in os.walk(d)])
    if results[d]==0: results[d]+=1
    total += results[d]

# Report sorted results:
for d in sorted(results, key=results.__getitem__, reverse=True):
    print('{:-<70s}{:->10d}'.format(d, results[d]))
print('TOTAL FILES: {:d}'.format(total))
