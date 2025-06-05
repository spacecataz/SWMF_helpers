#!/usr/bin/env python3

'''
Given a directory to ...

The result is written to JSON-headed ASCII; the default file name is
'iono_fluxes.txt'. The file can be read using spacepy as follows:

```
from spacepy import datamodel
data = datamodel.readJSONheadedASCII('iono_fluxes.txt')
```

'''

from spacepy.pybats import rim
import numpy as np
import os
from spacepy.datamodel import dmarray, SpaceData
import datetime
from argparse import ArgumentParser, RawDescriptionHelpFormatter

parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('folder', type=str, help='Directory to iono files')
parser.add_argument("-o", "--outfile", default='iono_fluxes.txt', help="Set " +
                    "output file name.  Defaults to 'iono_fluxes.txt'")

# Handle arguments:
args = parser.parse_args()


# Midpoint integration
def flux(var, theta, rads):
    '''
    Integrates the given flux 'var' over area using the midpoint method.
    It uses its ____ 'theta' and distance between measurements
    'rads' in the form [dTheta, dPhi]
    '''
    answer = var[i]*np.sin(theta)*np.square(6371000.0)*rads[0]*rads[1]
    return answer


# calculate netflux in gigawatts
def netflux(flux):
    '''
    Sums the integrated flux for the entire area and rounds the result to
    Gigawatts
    '''
    answer = np.sum(flux) / np.power(10, 6)
    return answer


# Setting directory to folder of run and getting '.gz' files to run
directory = args.folder
unfiltered_lst = sorted(os.listdir(directory))
lst = []
for L in unfiltered_lst:
    if os.path.splitext(L)[1] == '.gz':
        lst.append(L)
numfiles = len(lst)

# variables to calculate excluding 'time' and 'net_flux'
varis = ['diffuse', 'mono', 'elec', 'diffi', 'bbnd']

# initialize data
data = {}
netdata = SpaceData()
netdata['time'] = dmarray(np.zeros(numfiles, dtype=datetime.datetime))
for v in varis:
    netdata[v] = dmarray(np.zeros(numfiles))
netdata['net_flux'] = dmarray(np.zeros(numfiles))

# Analyze every file in directory
for L in range(numfiles):

    # read data in current file
    Iono = rim.Iono(f"{directory}/{lst[L]}")

    # initialize constants
    shape = (Iono.attrs['ntheta'], Iono.attrs['nphi'])
    length = shape[0] * shape[1]
    param = [np.radians(0.5), np.radians(1)]  # [dTheta, dPhi]
    theta = np.radians(Iono['n_theta'].flatten())

    # Initialize data vars for integration
    for k in varis:
        data[k] = dmarray(np.zeros(length))
        if k in ['diffi', 'bbnd', 'elec']:
            if k == 'elec' and 'n_e-flux-mono' in Iono:
                Kay = "n_e-flux-mono"
            else:
                Kay = f"n_e-flux-{k}"

            var = Iono[Kay].flatten()
            for i in range(length):

                data[k][i] = flux(var, theta[i], param)

    # Handle change from mono -> elec in simulation
    if 'n_e-flux-mono' in Iono:
        elec = Iono["n_e-flux-mono"].flatten()
    elif 'n_e-flux-elec' in Iono:
        elec = Iono["n_e-flux-elec"].flatten()

    diffe = Iono["n_e-flux-diffe"].flatten()

    # Calculate mono and diffuse portions of elec precipitation
    for i in range(length):

        if elec[i] >= 2*diffe[i]:
            data['mono'][i] = flux(elec, theta[i], param)
        elif elec[i] < 2*diffe[i]:
            data['diffuse'][i] = flux(elec, theta[i], param)
    for k in varis:
        netdata[k][L] = netflux(data[k])
    netdata['time'][L] = Iono.attrs['time']
    netdata['net_flux'][L] = np.sum([netdata['diffi'][L], netdata['elec'][L],
                                    netdata['bbnd'][L]])

netdata.toJSONheadedASCII(args.outfile)
