#!/usr/bin/env python3
'''
If FFMpeg is installed, it is simple to make movies out of a series
of PyBats PNG files.  Just use this script!  It's not just for PyBats;
it should work with any image series.

This script will combine any files you give it (usual Linux wildcard 
characters accepted) with the FFMpeg command to create a movie.

As an intermediate step to making the movie, this script places
the image files into /tmp/.  If this location is unavailable, use
the -t option to change it.

Usage:
make_movie.py [options] [files]

Options:

-h or -help:
     Print this help information.

-t=[directory]:
     Change the temporary directory.  Default is /tmp/

-r=[integer]:
     Set the frame rate for the input PNGs.  Default is 25.  The lower this
     number is, the slower the movie will progress through the PNGs.

-o=[file name]:
     Set the output file name.  Note that the extension sets the
     output format, so always include it!  .mpg, .mp4, and .avi
     all work well with the latter two being necessary for non-
     standard frame rates and the first two being preferrable for
     portability.  Default name is 'movie.mpg'.

'''

# This is not used right now:
#'''
#-b=[integer]:
#     Set the bit rate for ffmpeg.  Default is 2400 and is medium/low quality.
#
#'''

import sys
import glob
import os
import shutil

# Default/initial values:
files = []
rate=25
bps =2400
tmp = '/tmp/'
outfile = 'movie.mpg'

for option in sys.argv[1:]:
    # Handle options:
    if option[0] == '-':
        if option == '-h' or option == '-help':
            print(__doc__)
            exit()
        if option[0:2] == '-r':
            loc = option.rfind('=')
            rate = int(option[loc+1:])
        if option[0:2] == '-o':
            loc = option.rfind('=')
            outfile = option[loc+1:]
        if option[0:2] == '-t':
            loc = option.rfind('=')
            tmp = option[loc+1:]
            if tmp[-1] != '/':
                tmp = tmp + '/'
        if option[:2] == '-b':
            bps=int(option[option.rfind('=')+1:])
    # Search for files that match option.
    else:
        files = files + glob.glob(option)

if len(files) == 0:
    exit()

# Move files
for i, ifile in enumerate(files):
    shutil.copyfile(ifile, '%simg_%08d.png' % (tmp,i))

# Make movie
cmd='ffmpeg -r {} -i '.format(rate) + \
    '{}img_%08d.png -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -r 30 -pix_fmt yuv420p {}'.format(
        tmp, outfile)

print(cmd)

os.system(cmd)
#os.system('ffmpeg -s hd720 -r %i -i %simg_%%08d.png %s'
#          % (rate, tmp, outfile))

# Remove temp files.
for ifile in glob.iglob('/tmp/img_*.png'):
    os.remove(ifile)
