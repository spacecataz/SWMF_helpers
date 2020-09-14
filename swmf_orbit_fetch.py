#!/usr/bin/env python3

'''
Fetch orbit data from CDAweb for modern NASA missions and convert into
SWMF input files.
'''

from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
# Add arguments:
parser.add_argument('time', help='Start and end date to begin fetching in '+
                    'YYYYMMDD format.', type=str, nargs='*')
#parser.add_argument('start', help='Start date to begin fetching in YYYYMMDD'+
#                    ' format.', type=str)
#parser.add_argument('stop', help='Stop date of fetching in YYYYMMDD'+
#                    ' format.', type=str)
parser.add_argument('--save', '-s', help='Save intermediary files' +
                    ' Default behavior is to delete downloaded files.',
                    action='store_true')
parser.add_argument('--debug', '-d', help='Turn on debug output.',
                    action='store_true')
parser.add_argument('--verbose', '-v', help='Turn on verbose output',
                    action='store_true')
parser.add_argument('--plot', '-p', help='Generate plot of orbits',
                    action='store_true')
parser.add_argument('--info', '-i', help='Print available satellites and exit.',
                    action='store_true')
parser.add_argument('--sats', type=str, default='all',
                    help='Select what satellite orbits to fetch as a comma-'+\
                    'separated list.  Default is all satellites.')
parser.add_argument('--coord', '-c', help='Set output coordinate system. '+\
                    'Options are GSM, GSE, GEO, or SM.  Default is GSM.',
                    type=str, default='GSM')
parser.add_argument('--dt', type=float, default=15.0,
                    help='Set output frequency in SWMF #SATELLITE command.  '+
                    'Default value is 15s output.')

# Get args from caller, collect arguments into a convenient object:
args = parser.parse_args()

### MAIN IMPORTS HERE.
# Some take a while to load, we want no delay on printing help statements.
import sys, os, urllib
import datetime as dt
import xml.etree.ElementTree as ET

from spacepy.pycdf import CDF
from spacepy.pybats import SatOrbit
from spacepy.time import Ticktock
from spacepy.coordinates import Coords

#### Function definitions:
def plot_orbits(orbits, lim=[-20, 20], coord='GSM'):
    '''
    Given a list of satellite orbit objects, plot the orbit in a simple
    magnetosphere.
    '''

    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    
    from spacepy.empiricals import getMagnetopause
    from spacepy.pybats import add_planet
    from spacepy.pybats.dipole import b_line
    from spacepy.plot import style

    style()

    # Some parameters:
    ncol=8
    top =.88 if len(orbits)<(2*ncol+1) else .855
    
    # Set up plot:
    fig = plt.figure(figsize=(11,5))
    fig.subplots_adjust(top=top,bottom=0.11,left=0.09,
                        right=0.95,hspace=0.2,wspace=0.265)
    a1, a2 = fig.subplots(1,2)

    for a in (a1, a2):
        add_planet(a)
        a.set_aspect('equal')
        a.set_xlim(lim)
        a.set_ylim(lim)
        a.set_xlabel(f'{coord} X ($R_E$)', size=16)
        a.xaxis.set_major_locator(MultipleLocator(10))
        a.yaxis.set_major_locator(MultipleLocator(10))
    a1.set_ylabel(f'{coord} Y ($R_E$)', size=16)
    a2.set_ylabel(f'{coord} Z ($R_E$)', size=16)
        
    # equatorial plane- add magpause:
    pause = getMagnetopause( {'P':[1.5], 'Bz':[-2]})
    a1.plot(pause[0,:,0], pause[0,:,1], 'k--', alpha=.75, zorder=1)

    # meridional plane- add field lines:
    for L in [3, 6.6]:
        x,y,z = b_line(L,0,0, 101)
        a2.plot(x,z, 'k-', alpha=.5, zorder=1)
        x,y,z = b_line(-L,0,0, 101)
        a2.plot(x,z, 'k-', alpha=.5, zorder=1)

    # Add orbits:
    lines, labels = [], []
    for s in orbits:
        # Add lines:
        l=a1.plot(s['xyz'][0,:], s['xyz'][1,:])[0]
        l=a2.plot(s['xyz'][0,:], s['xyz'][2,:])[0]
        lines.append(l)
        labels.append(s.attrs['file'][:-4])
        
        # Annotate:
        for i in (0, -1):
            a1.plot(s['xyz'][0,i], s['xyz'][1,i], 'o', c=l.get_c())
            a1.annotate(f"{s['time'][i]:%Y-%m-%dT%H:%M}",
                        [s['xyz'][0,i],s['xyz'][1,i]], c=l.get_c(), size=8)
            a2.plot(s['xyz'][0,i], s['xyz'][2,i], 'o', c=l.get_c())
            a2.annotate(f"{s['time'][i]:%Y-%m-%dT%H:%M}",
                        [s['xyz'][0,i],s['xyz'][2,i]], c=l.get_c(), size=8)

    # Add legend:
    fig.legend(lines, labels, loc='upper left', fontsize=10, ncol=8)
    
    return fig
        
def fetch_cdf(sat, tstart, tstop, verbose=False):
    '''
    Create URL and fetch CDF file from CDAWeb's RESTful service.
    
    Returns the saved file name on disk on success, False on failure.
    '''

    # Constants:
    base = 'https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/'
    ns   = r'{http://cdaweb.gsfc.nasa.gov/schema}'
    
    url =  \
        f"{base}{sat['set']}/data/{tstart:%Y%m%dT%H%M%S}Z,"+\
        f"{tstop:%Y%m%dT%H%M%S}Z/{sat['var']}?format=cdf"

    if verbose:
        print(f'\tURL = {url}')

    # Request CDF and fetch XML response.
    with urllib.request.urlopen(url) as response:
        xml = ET.parse(response).getroot()

    # Search XML response for errors.
    if xml.findall(ns+'Error'):
        print('\tNo data for this satellite (error returned)')
        if verbose:
            for x in xml.findall(ns+'Error'):
                print(f'\t{x.text}')
        return False

    # Get URL of CDF file, download, save to disk.
    cdfurl = xml.find(ns+'FileDescription').find(ns+'Name').text
    filename = cdfurl.split('/')[-1]
    if verbose: print(f'\tCDF URL = {cdfurl}')
    urllib.request.urlretrieve(cdfurl, filename)

    # Return name of file to caller:
    return filename
            
#### MAIN SCRIPT FUNCTIONALITY:
    
# If we are in debug mode, we want lots and lots of output!
# Therefore, we turn on verbose mode and save mode, too:
if args.debug:
    args.verbose=True
    args.save=True

# Start by creating a dictionary of known sats and all information required
# to fetch and parse their orbit data from CDAWeb.
sats = {
    'cluster1':{'set':'C1_CP_FGM_SPIN', 'var':'sc_pos_xyz_gse__C1_CP_FGM_SPIN',
                'time':'time_tags__C1_CP_FGM_SPIN','coord':'GSE','units':'km'},
    'cluster2':{'set':'C2_CP_FGM_SPIN', 'var':'sc_pos_xyz_gse__C2_CP_FGM_SPIN',
                'time':'time_tags__C2_CP_FGM_SPIN','coord':'GSE','units':'km'},
    'cluster3':{'set':'C3_CP_FGM_SPIN', 'var':'sc_pos_xyz_gse__C3_CP_FGM_SPIN',
                'time':'time_tags__C3_CP_FGM_SPIN','coord':'GSE','units':'km'},
    'cluster4':{'set':'C4_CP_FGM_SPIN', 'var':'sc_pos_xyz_gse__C4_CP_FGM_SPIN',
                'time':'time_tags__C4_CP_FGM_SPIN','coord':'GSE','units':'km'},
    'goes06':{'set':'G6_K0_EPS', 'var':'SC_pos_sm',
                'time':'Epoch','coord':'GSM','units':'km'},
    'goes07':{'set':'G7_K0_EPS', 'var':'SC_pos_sm',
              'time':'Epoch','coord':'GSM','units':'km'},
    'goes08':{'set':'G8_K0_EP8', 'var':'SC_pos_sm',
              'time':'Epoch','coord':'GSM','units':'km'},
    'goes09':{'set':'G9_K0_EP8', 'var':'SC_pos_sm',
              'time':'Epoch','coord':'GSM','units':'km'},
    'goes10':{'set':'GOES10_EPHEMERIS_SSC', 'var':'XYZ_GSM',
              'time':'Epoch','coord':'GSM','units':'re'},
    'goes11':{'set':'GOES11_EPHEMERIS_SSC', 'var':'XYZ_GSM',
              'time':'Epoch','coord':'GSM','units':'re'},
    'goes12':{'set':'GOES12_EPHEMERIS_SSC', 'var':'XYZ_GSM',
              'time':'Epoch','coord':'GSM','units':'re'},
    'goes13':{'set':'GOES13_EPHEMERIS_SSC', 'var':'XYZ_GSM',
              'time':'Epoch','coord':'GSM','units':'re'},
    'goes14':{'set':'GOES14_EPHEMERIS_SSC', 'var':'XYZ_GSM',
              'time':'Epoch','coord':'GSM','units':'re'},
    'goes15':{'set':'GOES15_EPHEMERIS_SSC', 'var':'XYZ_GSM',
              'time':'Epoch','coord':'GSM','units':'re'},
    'goes16':{'set':'GOES16_EPHEMERIS_SSC', 'var':'XYZ_GSM',
              'time':'Epoch','coord':'GSM','units':'re'},
    'goes17':{'set':'GOES17_EPHEMERIS_SSC', 'var':'XYZ_GSM',
              'time':'Epoch','coord':'GSM','units':'re'},
    'geotail':{'set':'GE_OR_DEF', 'var':'GSM_POS',
              'time':'Epoch','coord':'GSM','units':'km'},
    'themisa':{'set':'THA_L1_STATE', 'var':'tha_pos_gsm',
               'time':'tha_state_epoch','coord':'GSM','units':'km'},
    'themisb':{'set':'THB_L1_STATE', 'var':'thb_pos_gsm',
               'time':'thb_state_epoch','coord':'GSM','units':'km'},
    'themisc':{'set':'THC_L1_STATE', 'var':'thc_pos_gsm',
               'time':'thc_state_epoch','coord':'GSM','units':'km'},
    'themisd':{'set':'THD_L1_STATE', 'var':'thd_pos_gsm',
               'time':'thd_state_epoch','coord':'GSM','units':'km'},
    'themise':{'set':'THE_L1_STATE', 'var':'the_pos_gsm',
               'time':'the_state_epoch','coord':'GSM','units':'km'},
    'rbspa':{'set':'RBSP-A-RBSPICE_LEV-3_ESRHELT', 'var':'Position_GSM',
             'time':'Epoch','coord':'GSM','units':'re'},
    'rbspb':{'set':'RBSP-B-RBSPICE_LEV-3_ESRHELT', 'var':'Position_GSM',
             'time':'Epoch','coord':'GSM','units':'re'},
}

# If info requested, print and close:
if args.info:
    print('Satellites that can be fetched with this script:')
    for s in sats: print(f'\t{s}')
    print('See documentation for info on expanding available sats.')
    sys.exit()

# Handle time arguments.  Note that we used nargs='*' so that users
# can print sat info w/o throwing errors.  This means that we need
# to check our arguments a bit better...
if len(args.time)!=2:
    print(__doc__)
    raise ValueError('Two time inputs required: '+
                     'start and stop time.  See documentation.')
try:
    # Try to parse date into datetime object.
    tstart = dt.datetime.strptime(args.time[0], '%Y%m%d')
    tstop  = dt.datetime.strptime(args.time[1],  '%Y%m%d')
except ValueError:  # Specify the type of exception to be specific!
    # If we can't, stop the program and print help.
    print('ERROR: Could not parse date!')
    print(__doc__)
    sys.exit()

if args.verbose: print(f'Searching from {tstart} to {tstop}...')
    
# Set what satellites to fetch:
sats_now = sats.keys() if args.sats=='all' else args.sats.split(',')
if args.debug: print('Working on sats:', sats_now)

# Build orbits for all satellites:
orbits=[]
for s in sats_now:
    print(f'Working on satellite {s}...')
    
    # Check for data, fetch CDFs:
    if args.verbose: print('\tObtaining orbit...')
    try:
        filename = fetch_cdf(sats[s], tstart, tstop, verbose=args.verbose)
    except urllib.error.HTTPError:
        print('\tNO DATA FOR THIS SATELLITE/DATE')
        continue

    # Open CDF file:
    if args.verbose: print(f'\tLoading {filename}...')
    cdf = CDF(filename)

    # Some debug information (useful when adding new sats):
    if args.debug:
        print('\tCDF Variables:')
        print(cdf)
    
    # Perform coordinate conversion:
    if args.coord == sats[s]['coord']:
        xyz = cdf[sats[s]['var']][...]
    else:
        if args.verbose:
            print(f"\tRotating from {sats[s]['coord']} to {args.coord}...")
        xyz_file = Coords(cdf[sats[s]['var']][...], sats[s]['coord'], 'car',
                          ticks = Ticktock(cdf[sats[s]['time']][...]) )
        xyz = xyz_file.convert(args.coord, 'car').data
                     
    # Convert to satellite orbit object:
    if args.verbose: print(f'\tConverting to SWMF input file...')
    out = SatOrbit()
    out['time'] = cdf[sats[s]['time']][...]
    out['xyz']  = xyz.transpose()

    # Convert units:
    if sats[s]['units']=='km':
        out['xyz']/=6371
        
    # Add attributes to file:
    out.attrs['file'] = f'{s}.sat'
    out.attrs['coor'] = sats[s]['coord']
    if out.attrs['coor'] == 'SM': out.attrs['coor']='SMG'
    tnow=dt.datetime.now()
    out.attrs['head'] =  [f'Created with swmf_orbit_fetch.py on '+
                          f'{tnow:%Y-%m-%d %H:%M%S}']

    if args.verbose: print(f'\tSaving to disk...')
    out.write()

    # Save orbit for plotting:
    orbits.append(out)

    # Remove CDF files
    if not args.save:
        if args.verbose: print(f'\tRemoving {filename}...')
        os.remove(filename)

# Write out the SWMF command to include satellites:
with open('sats.include','w') as f:
    f.write(f'#SATELLITE\n{len(orbits)}\t\t\tnSatellite\n')
    for orb in orbits:
        f.write('MHD RAY date time\tStringSatellite\n')
        f.write('-1\t\t\tDnOutput\n')
        f.write(f'{args.dt:.1f}\t\t\tDtOutput\n')
        f.write(f'{orb.attrs["file"]}\n')
        
print(f'Finished. Found {len(orbits)} of {len(sats_now)} requested spacecraft.')

if args.plot and bool(orbits):
    if args.verbose: print('\tCreating figure...')
    fig = plot_orbits(orbits)
    fig.savefig(f'./orbits_{tstart:%Y%m%d_%H%M%S}.png')
