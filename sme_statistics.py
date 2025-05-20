#!/usr/bin/env python3

'''
Given a txt file from the 'process_sme' script, calculate mean,
median, max, min, and perhaps standard deviation, correlation, RMSE, skew,
and kurtosis for a given time interval.

doc-string
hour window, calculate mean, median, max, min, and perhaps standard deviation,
correlation, RMSE, skew, and kurtosis
'''

# import re
# import datetime as dt
# import warnings

from matplotlib.dates import num2date  # date2num
from scipy import stats  # interpolate
import numpy as np
from spacepy.datamodel import SpaceData  # dmarray
import spacepy.datamodel as dm
from argparse import ArgumentParser, RawDescriptionHelpFormatter

parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('txtfile', type=str, help='Name of the data file ' +
                    'to calculate the statistics of the data.')
# parser.add_argument("-o", "--outfile", default='sm_indexes.txt', help="Set "
#                    "output file name.  Defaults to 'sm_indexes.txt'")

# Handle arguments:
args = parser.parse_args()
data = dm.readJSONheadedASCII(args.txtfile)

data['time'] = num2date(data['time'])

vars = ['SML', 'SWMFL', 'mlat_L', 'SMLmlat']  # ... for testing
for v in vars:
    data[v] = SpaceData(Values=data[v])
    Vals = 'Values'
    data[v] = SpaceData(data[v], Max=max(data[v][Vals]))
    data[v] = SpaceData(data[v], Min=min(data[v][Vals]))
    data[v] = SpaceData(data[v], Mean=(np.mean(data[v][Vals]) * 1.))
    data[v] = SpaceData(data[v], Median=np.median(data[v][Vals]))
    data[v] = SpaceData(data[v], Standard_Deviation=np.std(data[v][Vals]))
    data[v] = SpaceData(data[v], Skew=stats.skew(data[v][Vals]))
    data[v] = SpaceData(data[v], Kurtosis=stats.kurtosis(data[v][Vals]))
    data[v] = SpaceData(data[v], Min=min(data[v][Vals]))
