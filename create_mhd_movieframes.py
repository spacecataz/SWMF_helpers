#!/usr/bin/env python3
'''
This script generates a series of PNG figures characterizing a SWMF simulation.
Each plot is arranged as follows:
       ------------------------------------------------------------
       |  ----------------------        -----------------------   |
       |  |                    |        |                     |   |
       |  |        Y=0         |        |        Z=0          |   |
       |  |      PRESSURE      |        |      DENSITY        |   |
       |  |                    |        |                     |   |
       |  |                    |        |                     |   |
       |  |                    |        |                     |   |
       |  ----------------------        -----------------------   |
       | -------------------------------------------------------  |
       | |                                                     |  |
       | |              DST RESULTS FOR EVENT                  |  |
       | |                                                     |  |
       | -------------------------------------------------------  |
       | -------------------------------------------------------  |
       | |                                                     |  |
       | |                IMF DATA FOR EVENT                   |  |
       | |                                                     |  |
       | -------------------------------------------------------  |
       ------------------------------------------------------------

A vertical line is placed in the IMF data indicating the epoch of each plot.
Input files are searched for in PWD only.  Results are placed in ./moviepngs
unless otherwise indicated.

This script should be run from the GM directory of an SWMF results directory
(i.e., one created by a call to the Postproc.pl script).
'''

import os
from sys import stdout
from glob import glob
import datetime as dt
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import multiprocessing as mp

import numpy as np
import matplotlib
from dateutil.parser import parse
import matplotlib.pyplot as plt

import spacepy.pybats as pb
import spacepy.pybats.bats as pbs
from spacepy.plot import applySmartTimeTicks, style

# Required defaults:
start = False
end = False

# Set up argruments:
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-o", "--outdir", default='mhd_pngs', help="Set output " +
                    "directory for plots.  Defaults to ./mhd_pngs")
parser.add_argument("-m", "--mag", action="store_true",
                    help="Add magnetic field lines to y=0 slice")
parser.add_argument("-obs", "--obsdst", help="Fetch and plot observed dst on" +
                    " top of modeled values", action="store_true")
parser.add_argument("--debug", help="Only the first plot is created.  " +
                    "It is shown to the user and not saved.  Additional " +
                    "information is printed to screen.", action="store_true")
parser.add_argument("-n", "--nthread", type=int, default=1,
                    help="If >1, plots will be made in parallel using N " +
                    "threads")
parser.add_argument("-t", "--title", default='', help="Plot title, in quotes")
parser.add_argument("-s", "--sats", action="store_true", help="Add virtual " +
                    "satellite positions and names to plot using all " +
                    "sat*.sat files found")
parser.add_argument("-i", "--imf", default='', help="Specify the IMF file.  " +
                    "If none given, the script will attempt to find one in " +
                    "the pwd or up one directory.")
parser.add_argument("-by", "--by", action='store_true', help="Turn on " +
                    "plotting of IMF BY.  Default is to not show BY.")
parser.add_argument("-v1", "--var1", default='p', help="Set the variable to " +
                    "plot in the Y=0 frame.  Defaults to 'p' for pressure")
parser.add_argument("-v2", "--var2", default='rho', help="Set the variable " +
                    "to plot in the Z=0 frame.  Defaults to 'rho' for total" +
                    " mass density")
parser.add_argument("-c1", "--cmap1", default='plasma', help='Set the ' +
                    'Matplotlib color map for the Y=0 frame.  Any MPL table ' +
                    'name can be used; defaults to "plasma"')
parser.add_argument("-c2", "--cmap2", default='viridis', help='Set the ' +
                    'Matplotlib color map for the Z=0 frame.  Any MPL table ' +
                    'name can be used; defaults to "viridis"')
parser.add_argument('-start', default=None, help="Set the start date and " +
                    "time of the simulation using the format " +
                    "YYYY-MM-DDTHH:MN:SS. If not present, an attempt will be" +
                    " made to obtain it from local PARAM files.")
parser.add_argument('-end', default=None,
                    help="Same as -start but for end time")
parser.add_argument('-z1', '--zlim1', default=(.01, 200), help='Set the ' +
                    'color bar range for the Y=0 slice using two numbers, ' +
                    'space separated (e.g., -z1 .01 200).  Log scale is used!',
                    type=float, nargs=2)
parser.add_argument('-z2', '--zlim2', default=(.1, 200), help='Set the color' +
                    ' bar range for the Z=0 slice using two numbers, space ' +
                    'separated (e.g., -z2 .1 200).  Log scale is used!',
                    type=float, nargs=2)
parser.add_argument('-xrng', '--xrng', default=(-35, 15), help='Set the ' +
                    'spatial extent in the GSM X direction via two space-' +
                    'separated numbers with no spaces (e.g., -xrng -10 10)',
                    type=float, nargs=2)
parser.add_argument('-yrng', '--yrng', default=(-22.5, 22.5),
                    help='Set the spatial extent in the GSM Y direction via ' +
                    'two space-separated numbers with no spaces ' +
                    '(e.g., -yrng -10 10)', type=float, nargs=2)
parser.add_argument('-zrng', '--zrng', default=(-22.5, 22.5),
                    help='Set the spatial extent in the GSM Z direction via ' +
                    'two space-separated numbers with no spaces ' +
                    '(e.g., -zrng -10 10)', type=float, nargs=2)
parser.add_argument('--info', action="store_true", help="Get information " +
                    "about simulation & available plot variables, " +
                    "but do not plot.")

# Handle arguments:
args = parser.parse_args()

# Some convenience variables from our parsed arguments:
imffile = os.path.expanduser(args.imf)
v1, v2 = args.var1, args.var2
cmap1, cmap2 = args.cmap1, args.cmap2
start, end = args.start, args.end

if args.debug:
    print('\t<<<OPERATING IN DEBUG MODE!>>>')
else:
    matplotlib.rcParams['figure.max_open_warning'] = 1E6
    matplotlib.use('Agg')

# Update to Spacepy Style
style()


def get_start_time(prefix):
    '''
    Try to find start time in PARAM files.  Use prefix to determine
    where PARAMs should be.
    '''

    from spacepy.pybats import parse_filename_time

    starttime = False

    # Start by searching for PARAM.in files in PWD and, if necessary,
    # up one or two directories.
    params = glob('PARAM*')
    if prefix == './GM/IO2/':
        params += glob('../../PARAM*')
    if prefix == './GM/':
        params += glob('../PARAM*')

    # Flip through files searching for #STARTTIME
    for p in params:
        f = open(p, 'r')
        line = f.readline()
        while line:
            if 'STARTTIME' in line:
                starttime = dt.datetime(
                    int(f.readline()[:4]), int(f.readline()[:2]),
                    int(f.readline()[:2]), int(f.readline()[:2]),
                    int(f.readline()[:2]), int(f.readline()[:2]))
                return starttime
            line = f.readline()

    # Next, attempt to use the file names themselves:
    y_file = glob(prefix+'/y??_mhd*.out')[0]
    i_iter, runtime, starttime = parse_filename_time(y_file)

    return starttime


# Other constants of interest:
bzcolor = '#3333CC'
bycolor = '#9999FF'
pdcolor = '#CC3300'

# Convert ranges from text to numbers:
prng = np.array(args.zlim1)
drng = np.array(args.zlim2)

xrng = np.array(args.xrng)
yrng = np.array(args.yrng)
zrng = np.array(args.zrng)

if args.debug:
    print("Scales and ranges:")
    print("\tZ1   = [{0[0]},{0[1]}]".format(prng))
    print("\tZ2   = [{0[0]},{0[1]}]".format(drng))
    print("\tXrng = [{0[0]},{0[1]}]".format(xrng))
    print("\tYrng = [{0[0]},{0[1]}]".format(yrng))
    print("\tZrng = [{0[0]},{0[1]}]".format(zrng))


# Create a folder for holding our PNGS.
if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)

# Determine where to look for GM output files.
prefix = './'
if os.path.exists('./GM/IO2'):
    prefix += 'GM/IO2/'
elif os.path.exists('./GM'):
    prefix += 'GM/'
print('\tSearching in folder {}'.format(prefix))

# Set start time.  Search for it as necessary.
if start:
    # Convert to datetime.
    start = parse(start)
else:
    # Try to find in param files.
    start = get_start_time(prefix)
# If we fail to find start time, stop code.
if not start:
    raise(ValueError('Simulation start time not found, use -start option.'))
print('\tStart time: {}'.format(start.isoformat()))

# No imf file given?  Let's try to find one.
if not imffile:
    if glob('imf*.dat'):
        imffile = glob('imf*.dat')[0]
    elif glob('IMF*.dat'):
        imffile = glob('IMF*.dat')[0]
    elif glob('../imf*.dat'):
        imffile = glob('../imf*.dat')[0]
    elif glob('../IMF*.dat'):
        imffile = glob('../IMF*.dat')[0]
    else:
        raise(ValueError("Could not find IMF file.  Use -i flag."))
    print('\tIMF File: using {}'.format(imffile))

# Obtain imf and log files.  Skip steady-state log files.
imf = pb.ImfInput(imffile)
for logname in glob(prefix+'log_*.log'):
    log = pbs.BatsLog(logname, starttime=start)
    if log['runtime'][-1] != log['runtime'][0]:
        break

# Did we find a log file?
if 'log' not in locals():
    raise(ValueError("Did not find MHD Logfile (required)."))
print('\tLog File: using {}'.format(logname))

# Get observed DST.  Raise exception if we cannot.
if args.obsdst:
    if not log.fetch_obs_dst():
        raise ValueError('Failed to obtain observed Dst.  Is KyotoWDC down?')

imf.calc_pram()

if end:
    end = parse(end)
else:
    end = log['time'][-1]
trange = [start, end]

sats = []
if args.sats:
    for s in glob(prefix+"*.sat"):
        sats.append(pbs.VirtSat(s))
if sats:
    print('\tVirtSats: found {} virtual satellites'.format(len(sats)))

if args.nthread > 1:
    print('\tParallel: using {} threads.'.format(args.nthread))


# Create function for generating plots.
def plot_results(fileY, fileZ, iFile=0, nFiles=1):
    '''
    Create a nice plot of WTF is happening.
    Kwargs iFile and nFiles are for reporting progress.
    '''

    # Get start times from file names:
    time_y = pb.parse_filename_time(fileY)
    time_z = pb.parse_filename_time(fileZ)

    # Check that times match:
    for t1, t2 in zip(time_y, time_z):
        if t1 != t2:
            print(f'File times do not match!\n\tY-File: {t1}\n\tZ-File: {t2}')
            raise ValueError('File time mismatch')

    # Use file time info get actual time:
    if time_y[-1]:
        tNow = time_y[-1]      # Datetime in file name.
    elif time_y[1]:
        t_delta = dt.timedelta(seconds=time_y[1])
        tNow = start + t_delta  # Runtime in file name.
    else:
        tNow = start           # Only iteration given.

    # Check if files are outside time range.
    if tNow < trange[0]:
        return 0
    if tNow > trange[-1]:
        return 0

    # Open files & calculate important values:
    mY = pbs.Bats2d(fileY)
    mZ = pbs.Bats2d(fileZ)
    mY.calc_ndens()
    mZ.calc_ndens()
    mY.calc_utotal()
    mZ.calc_utotal()

    if mZ.attrs['runtime'] != mY.attrs['runtime']:
        print('{}\n{}\n'.format(fileY, fileZ))
        raise(ValueError('Times of files do not line up!'))

    fig = plt.figure(figsize=[11, 8.5])
    fig.subplots_adjust(left=0.08, right=.96, bottom=0.06, top=.93)

    a1 = fig.add_subplot(221)
    stuff = mY.add_contour('x', 'z', v1, 91, dolog=True, xlim=xrng, ylim=zrng,
                           add_cbar=True, target=a1, zlim=prng, cmap=cmap1)
    a1.set_title(f'{pb.mhdname_to_tex(v1)} ({mY[v1].attrs["units"]})')

    stuff[-1].set_label('')
    if args.mag:
        mY.add_b_magsphere(a1, DoOpen=True)
    a2 = fig.add_subplot(222)
    stuff = mZ.add_contour('x', 'y', v2, 91, dolog=True, xlim=xrng, ylim=yrng,
                           add_cbar=True, target=a2, zlim=drng, cmap=cmap2)
    a2.set_title(f"{pb.mhdname_to_tex(v2)} ({mZ[v2].attrs['units']})")
    stuff[-1].set_label('')

    for s in sats:
        s.add_sat_loc(tNow, a1, plane='XZ', dolabel=True, color='r', size=9)
        s.add_sat_loc(tNow, a2, plane='XY', dolabel=True, color='r', size=9)

    a3, a4 = fig.add_subplot(413), fig.add_subplot(414)
    a5 = a4.twinx()
    pos3, pos4 = a3.get_position(), a4.get_position()
    a3.set_position([pos3.xmin, pos3.ymin, pos3.width*.95, pos3.height])
    a4.set_position([pos4.xmin, pos4.ymin, pos4.width*.95, pos4.height])
    a5.set_position([pos4.xmin, pos4.ymin, pos4.width*.95, pos4.height])

    log.add_dst_quicklook(target=a3, plot_obs=args.obsdst,
                          obs_kwargs={'c': 'k', 'ls': '-'}, color='r')
    a3.legend(loc='best', fontsize=10)
    a3.set_xlabel('')
    applySmartTimeTicks(a3, trange, dolabel=False)
    a3.grid(False, axis='y')
    a3.set_ylabel('D$_{ST}$ ($nT$)')
    a3.hlines(0, trange[0], trange[1], color='grey', linestyles='dashed')
    ymin, ymax = a3.get_ylim()
    a3.vlines(tNow, ymin, ymax, linestyles='solid', color='g', linewidths=2.)
    a3.set_ylim([ymin, ymax])

    a4.plot(imf['time'], imf['bz'],   color=bzcolor, label='B$_Z$')
    if args.by:
        a4.plot(imf['time'], imf['by'], color=bycolor, label='B$_Y$')
        a4.legend(loc='best', ncol=2)
    a5.plot(imf['time'], imf['pram'], color=pdcolor, lw=1.5)
    applySmartTimeTicks(a4, trange, dolabel=True)
    a4.grid(False, axis='y')
    a4.hlines(0, trange[0], trange[1], color=bzcolor, linestyles='dashed')
    a4.set_ylabel('IMF '+(not args.by)*'B$_{Z}$ ' + '($nT$)', color=bzcolor)
    a5.set_ylabel('P$_{dyn}$ ($nPa$)',  color=pdcolor)
    a5.grid(False, axis='y')
    ymin, ymax = a4.get_ylim()
    a4.vlines(tNow, ymin, ymax, linestyles='solid', color='g', linewidths=2.)
    a4.set_ylim([ymin, ymax])

    fig.suptitle(args.title, size=18)

    # Save the figure (normal operation) or show figure (debug mode).
    if args.debug:
        plt.show()
    else:
        fig.savefig(args.outdir+'/mhd{:%Y%m%d_%H%M%S}.png'.format(tNow))

    # Clear and close figure.
    if not args.debug:
        fig.clf()
        plt.close(fig)

    percent_done = float(iFile+1)/float(nFiles)*100.0
    return percent_done


def print_progress(percent):
    stdout.write(19*'\b')
    stdout.write('\t{:5.1f}% complete...'.format(percent))
    stdout.flush()


if __name__ == '__main__':
    # Find files, sort them.
    yFiles = glob(prefix+'y??_mhd*_[te]*.out')
    zFiles = glob(prefix+'z??_mhd*_[te]*.out')
    yFiles.sort()
    zFiles.sort()
    if not (yFiles and zFiles):
        raise(ValueError("Could not find any Y or Z slice files."))

    if args.debug:
        yFiles, zFiles = yFiles[:1], zFiles[:1]

    if args.info:
        mhd = pbs.Bats2d(yFiles[0])
        mhd.calc_ndens()
        mhd.calc_utotal()
        mhd.tree()
        exit()

    # pool.map(plot_results, zip(yFiles, zFiles)[:3])
    # Loop over files, making a plot for each set.
    nFiles = len(yFiles)
    print('\tFiles: working on {} files...'.format(nFiles))
    stdout.write('\t{:5.1f}% complete...'.format(0))

    if args.nthread > 1:
        # Create pool of workers.
        pool = mp.Pool(args.nthread)

        for i, files in enumerate(zip(yFiles, zFiles)):
            pool.apply_async(plot_results, args=(files[0], files[1]),
                             kwds={'iFile': i, 'nFiles': nFiles},
                             callback=print_progress)

        pool.close()
        pool.join()

    else:
        for i, files in enumerate(zip(yFiles, zFiles)):
            per = plot_results(files[0], files[1], i, nFiles)
            print_progress(per)

    print('\nALL DONE\n')
