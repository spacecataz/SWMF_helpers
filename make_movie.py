#!/usr/bin/env python3
'''
If FFMpeg is installed, it is simple to make movies out of a series
of PNG (or similar) files.  Just use this script!  It's not just for PyBats;
it should work with any image series.

This script will combine any files you give it (usual Linux wildcard
characters accepted) with the FFMpeg command to create a movie.

As an intermediate step to making the movie, this script places
the image files into /tmp/.  If this location is unavailable, use
the -t option to change it.

'''

import glob
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Potentially no longer used:
# import shutil

# Start by configuring the argparser:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("files", nargs='+', help="File glob pattern for " +
                    "ordered images to convert (e.g., images_*.png).")
parser.add_argument("-o", "--outfile", default='movie.mp4', help="Set the " +
                    "output file name.  Note that the extension sets the " +
                    "output format, so always include it!  .mpg, .mp4, and " +
                    ".avi all work well with the latter two being necessary " +
                    "for non-standard frame rates and the first two being " +
                    "preferrable for portability.")
parser.add_argument("-t", "--tempdir", default='/tmp/', help="Set the " +
                    "temporary directory to which files are copied.")
parser.add_argument("-r", "--rate", default=25, type=int, help="Set the " +
                    "frame rate for the input PNGs.  Default is 25.  The " +
                    "lower this number is, the slower the movie will " +
                    "progress through the PNGs.")
parser.add_argument("-n", "--nloop", default=0, type=int, help="Set number " +
                    "of times the video loops. 0 is no looping, -1 is " +
                    "infinite looping.")
parser.add_argument("--debug", default=False, action='store_true',
                    help="Turn on debugging mode.")
args = parser.parse_args()

# Default/initial values not handled above:
bps = 2400

if len(args.files) == 0:
    print("No files found to convert. Check input syntax.")
    exit()

# Create sorted list of file names to use as input:
with open('.make_movie_input.txt', 'w') as outfile:
    args.files.sort()
    for f in args.files:
        outfile.write(f"file '{f}'\n")

# Move files (no longer used?)
# for i, ifile in enumerate(args.files):
#     shutil.copyfile(ifile, '%simg_%08d.png' % (tmp,i))

# Make movie
cmd = 'ffmpeg -stream_loop {:d} -r {} '.format(args.nloop, args.rate) + \
      '-f concat -i .make_movie_input.txt -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ' + \
      '-c:v libx264 -pix_fmt yuv420p {}'.format(args.outfile)
# args.files -> args.tempdir

if args.debug:
    print(cmd)

# Execute command:
os.system(cmd)

# Remove temp file:
if not args.debug:
    os.remove('.make_movie_input.txt')
