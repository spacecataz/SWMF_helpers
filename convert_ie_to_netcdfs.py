#!/usr/bin/env python3

'''
This script converts IE/Ridley_serial output into a form that can be used by
the WAM-IPE team.  The hemispheres are combined into a single grid, and
coordinates converted from geomagnetic to geographic.  The new output is 
written to a NetCDF file.
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as np
from scipy.io.netcdf import netcdf_file
from spacepy.pybats import rim

def combine(ie):
    '''
    Given an ionospheric slice, rotate to geographic coordinates.
    '''

    from spacepy.datamodel import dmarray

    nHem = ie.attrs['ntheta']      # npoints in one hemi.
    nLat = 2*ie.attrs['ntheta']-1  # combine northern and southern hemi.
    nLon = ie.attrs['nphi']
    
    # New coords- 1D for geomagnetic (uniform), 2D for geographic (non-uniform):
    ie['lat'] = np.zeros(nLat)
    ie['lon'] = np.zeros(nLon)

    # We can immediately fill the 'lat' and 'lon' arrays:
    ie['lat'][:nHem] = 90-ie['n_theta'][:, 0]
    ie['lat'][nHem:] = 90-ie['s_theta'][1:,0]
    ie['lon'] = ie['n_psi'][0,:]
    
    # Combine the hemispheres of the x,y,z-directed fields in the new arrays.
    ie['phi'] = dmarray( np.zeros((nLat,nLon)), {'units': 'kV'})

    # combine hemispheres for old values:
    ie['phi'][:nHem,:] = ie['n_phi']
    ie['phi'][nHem:,:] = ie['s_phi'][1:,:]

def create_netcdf(iono, outdir=None):
    '''
    Given a rotated Ionosphere object, spit it out to NetCDF.
    '''

    combine(iono)
    
    # Create filenames for geomagnetic and geograpic results:
    outfile1 = 'potential_SMG_{:%Y%m%d_%H%M%S}.nc'.format(iono.attrs['time'])

    # Add output directory if given:
    if outdir:
        outfile1 = outdir + '/' + outfile1
    
    # Set up files & meta data:
    out1 = netcdf_file(outfile1, 'w')
    out1.time   = '{}'.format(iono.attrs['time'])
    out1.coords = 'Lat/Lon in SM coordinates.'

    # Create dimensions:
    out1.createDimension('lat', iono['lat'].size)
    out1.createDimension('lon', iono.attrs['nphi']  )

    # Add lat/lon to each:
    lat1 = out1.createVariable( 'lat', 'f', ['lat'] )
    lon1 = out1.createVariable( 'lon', 'f', ['lon'] )
    lat1.units='degress';  lon1.units='degrees'
    lat1[:] = iono['lat']; lon1[:] = iono['lon']

    # Add potential:
    e1 = out1.createVariable('potential', 'f', ('lat', 'lon') )
    e1[:,:]  = iono['phi']; e1.units = iono['phi'].attrs['units']
        
    out1.close()


def plot_rotated(ie):
    '''
    Create a series of plots to examine the rotation:
    '''

    import matplotlib.pyplot as plt
    
    for d in 'nev':

        # Get zlimit/plot range for each plot:
        z = 0
        for x in 'nev': z = max(z, np.abs(ie['e'+x]).max())
        z *= .75
        
        # Create figure:
        fig = plt.figure( figsize=(8.5,11) )
        
        # Add polar plots at top:
        out = ie.add_cont('n_e'+d, n=31, target=fig, loc=321, maxz=z, add_cbar=True)
        out = ie.add_cont('s_e'+d, n=31, target=fig, loc=322, maxz=z, add_cbar=True)

        # Add maps below:
        a3, a4 = fig.add_subplot(312), fig.add_subplot(313)
        levels = out[2].levels
        
        a3.contourf(ie['lon'], ie['lat'], np.array(ie['e'+d]),    levels=levels,
                    cmap=rim.get_iono_cb())
        a4.tricontourf(ie['glon'].flatten(), ie['glat'].flatten(),
                       np.array(ie['e'+d].flatten()), levels=levels,
                       cmap=rim.get_iono_cb())

        a3.grid(); a4.grid()
        
        a3.set_ylabel('SM Lat')
        a4.set_ylabel('Geo Lat')
        
        
def test_rotation():
    '''
    '''

    import matplotlib.pyplot as plt
    
    # An example file to work with:
    ie = rim.Iono('Wam_20150316/IE/it150317_142400_000.idl')

    # Rotate file:
    rotate(ie)

    # Save NetCDF:
    create_netcdf(ie)
    
    # Get maximum field:
    z = 0
    for x in 'nev':
        z = max(z, np.abs(ie['n_e'+x]).max())
    
    # Create a big plot:
    fig = plt.figure()
    ie.add_cont('n_phi', add_cbar=True, target=fig, loc=121)
    ie.add_cont('n_jr',  add_cbar=True, target=fig, loc=122)
    fig.suptitle(ie.attrs['time'])
    fig = plt.figure( figsize=(11,8) )
    ie.add_cont('n_ex', target=fig, loc=231, maxz=z, add_cbar=True)
    ie.add_cont('n_ey', target=fig, loc=232, maxz=z, add_cbar=True)
    ie.add_cont('n_ez', target=fig, loc=233, maxz=z, add_cbar=True)
    ie.add_cont('n_en', target=fig, loc=234, maxz=z, add_cbar=True)
    ie.add_cont('n_ee', target=fig, loc=235, maxz=z, add_cbar=True)
    ie.add_cont('n_ev', target=fig, loc=236, maxz=z, add_cbar=True)

    # Create moar plots:
    plot_rotated(ie)

    return ie

if __name__ == '__main__':


    # Start by setting up and parsing arguments.
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs='+',
                        help="IE files to convert.  Accepts Unix wildcards.")
    parser.add_argument("-o", "--outdir", default='./potential_netcdfs/',
                        help="Set output directory for NetCDF files. " +
                        "Defaults to ./potential_netcdfs")
    args = parser.parse_args()
    
    import os
    from spacepy.pybats.rim import Iono

    if not os.path.exists(args.outdir): os.mkdir(args.outdir)
    
    for ie in args.files:
        raw = Iono(ie)
        print(f"Working on {ie}...")
        create_netcdf(raw, outdir=args.outdir)
