#!/usr/bin/env python3

'''
This script unpacks SWMF-formatted ".outs" files into individual ".out"
components.

The SWMF can produce ".out" files, which represent output from a single point
in time. These can be concatenated into ".outs" files (note the "s") which
contain all data from all individual ".out" files. On certain occasions, it
is desireable to unpack ".outs" files into a series of ".out" files. This
script accomplishes that.

'''

import os
import sys
import re
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("files", nargs='+', help="*.outs files to unpack. " +
                    "Can be explicit files or a unix wildcard.")
parser.add_argument('--remove', '-rm', help="Remove original .outs file " +
                    'after concatenation leaving only the .outs.',
                    default=True, action='store_true')
parser.add_argument('--verbose', '-v', default=False, action='store_true',
                    help='Turn on verbose mode.')
parser.add_argument('--debug', '-d', default=False, action='store_true',
                    help='Turn on debug information.')

