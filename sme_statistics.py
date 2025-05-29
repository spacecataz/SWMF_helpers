#!/usr/bin/env python3

'''
Given a txt file and time interval from the 'process_sme' script, calculate
mean, median, max, min, and perhaps standard deviation, correlation, RMSE,
skew, and kurtosis for a given time interval.

Input format: data.txt -s %Y-%m-%dT%H:%M:%S -e %Y-%m-%dT%H:%M:%S
'''

from matplotlib.dates import date2num  # num2date
from scipy import stats  # interpolate
import numpy as np
import spacepy.datamodel as dm
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import datetime


parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('txtfile', type=str, help='Name of the data file ' +
                    'to calculate the statistics of the data.')
parser.add_argument("-s", "--starttime", default=None, help='Start time' +
                    ' using dt(year, month, day, hour, minute, second). ' +
                    'Default to beginning of file.')
parser.add_argument("-e", "--endtime", default=None, help='End time using ' +
                    'dt(year, month, day, hour, minute, second). ' +
                    'Defaults to end of file.')
# parser.add_argument("-o", "--outfile", default='sm_indexes.txt', help="Set "
#                    "output file name.  Defaults to 'sm_indexes.txt'")

# Handle arguments:
args = parser.parse_args()
data = dm.readJSONheadedASCII(args.txtfile)
if args.starttime is None:
    start_time = None
else:
    start_time = datetime.datetime.strptime(args.starttime, '%Y-%m-%dT%H:%M:%S'
                                            )

if args.endtime is None:
    end_time = None
else:
    end_time = datetime.datetime.strptime(args.endtime, '%Y-%m-%dT%H:%M:%S')
# fix dates #### NOTE: Should be adressed later
# data['time'] = num2date(data['time'])
# another thing for the time issue


# define the class
class SMinterval:
    swmfu = {}
    swmfl = {}
    '''
    Given a mag_grid data object, a start time, and and end time;
    returns statistics on the mlat and mlt data for both the smu and
    SML.

    If only one or no times are given, defaults to beginning or ending
    time respectively.

    SMinterval(mag_object, datetime.datetime(start_time),
               datetime.datetime(end_time))
    '''
    def __init__(self, mag_object, start=None, end=None):
        self.start = start
        self.end = end
        self.swmfl['mlat'] = self.calc_from_mags(mag_object, 'L', 'mlat')
        self.swmfl['mlt'] = self.calc_from_mags(mag_object, 'L', 'mlt')
        self.swmfu['mlat'] = self.calc_from_mags(mag_object, 'U', 'mlat')
        self.swmfu['mlt'] = self.calc_from_mags(mag_object, 'U', 'mlt')

    def calc_from_mags(self, mags, sm_index, value):
        # Turns given time into numerical values for rounding
        # If no time input is given defaults to beginning/ending time
        if self.start is None:
            start_time = mags['time'][0]
        else:
            start_time = date2num(self.start)

        if self.end is None:
            end_time = mags['time'][-1]
        else:
            end_time = date2num(self.end)

        # Using the times pull the correct data for analysis
        for i in range(len(mags['time'])):
            if round(mags['time'][i], 4) == round(start_time, 4):
                self.sframe = i
            elif round(mags['time'][i], 4) == round(end_time, 4):
                self.eframe = i + 1
            else:
                pass

        # Get correct range of data set fors analysis
        index = value + '_' + sm_index
        data = mags[index][self.sframe:self.eframe]
        list = ['SML', 'SMLmlat', 'SMLmlt', 'SMU', 'SMUmlat', 'SMUmlt',
                'SWMFL', 'SWMFU', 'lon_L', 'lon_U', 'mlat_L', 'mlat_U',
                'mlt_L', 'mlt_U', 'time']
        twags = {}
        for L in list:
            twags[L] = mags[L][self.sframe:self.eframe]

        # Calculate the statistics
        return SMdata(data, twags, index)


class SMdata:
    '''
    Given a dataset, returns statistics for analyis.

    Stats calculated:
    Max
    Min
    Mean
    Median
    Standard Deviation
    Skew
    Kurtosis
    Correlation Coefficient
    RSME
    '''
    def __init__(self, data_in, mags, index):
        self.data = data_in
        self.stats = {}
        self.correlation = {}
        self.calc_stats(mags, index)
        self.calc_corr(mags)

    def __setitem__(self, key, value):
        if isinstance(key, int) | isinstance(key, slice):
            self.data[key] = value
        elif isinstance(key, str):
            self.stats[key] = value

    def __getitem__(self, item):
        if isinstance(item, int) | isinstance(item, slice):
            return self.data[item]
        elif isinstance(item, str):
            return self.stats[item]

    def calc_stats(self, mags, index):
        # Calculate data statists
        self.stats['mean'] = np.mean(self.data)
        self.stats['median'] = np.median(self.data)
        self.stats['max'] = max(self.data)
        self.stats['min'] = min(self.data)
        self.stats['standard deviation'] = np.std(self.data)
        self.stats['skew'] = stats.skew(self.data)
        self.stats['kurtosis'] = stats.kurtosis(self.data)

        # Calculate RSME and Correlation Coefficient between
        # SWMF and SuperMAG values
        if index == 'mlat_L':
            MSE = np.square(np.subtract(mags['SMLmlat'], self.data)).mean()
            RMSE = np.sqrt(MSE)
            self.stats['correlation coefficient'] = np.corrcoef(self.data,
                                                                mags['SMLmlat']
                                                                )[0, 1]
        elif index == 'mlt_L':
            MSE = np.square(np.subtract(mags['SMLmlt'], self.data)).mean()
            RMSE = np.sqrt(MSE)
            self.stats['correlation coefficient'] = np.corrcoef(self.data,
                                                                mags['SMLmlt']
                                                                )[0, 1]
        elif index == 'mlat_U':
            MSE = np.square(np.subtract(mags['SMUmlat'], self.data)).mean()
            RMSE = np.sqrt(MSE)
            self.stats['correlation coefficient'] = np.corrcoef(self.data,
                                                                mags['SMUmlat']
                                                                )[0, 1]
        elif index == 'mlt_U':
            MSE = np.square(np.subtract(mags['SMUmlt'], self.data)).mean()
            RMSE = np.sqrt(MSE)
            self.stats['correlation coefficient'] = np.corrcoef(self.data,
                                                                mags['SMUmlt']
                                                                )[0, 1]
        self.stats['RSME'] = RMSE

    def calc_corr(self, mags):
        # Correlation Coefficient compared to every value
        list = ['SML', 'SMLmlat', 'SMLmlt', 'SMU', 'SMUmlat', 'SMUmlt',
                'SWMFL', 'SWMFU', 'lon_L', 'lon_U', 'mlat_L', 'mlat_U',
                'mlt_L', 'mlt_U', 'time']
        for L in list:
            self.correlation[L] = np.corrcoef(self.data, mags[L])


# Output data
Stats = SMinterval(data, start=start_time, end=end_time)

# Write stats to new file
statistics = ['mean', 'median', 'max', 'min', 'standard deviation', 'skew',
              'kurtosis', 'correlation coefficient', 'RSME']
with open('Stats_' + args.txtfile, 'w') as f:
    f.write(datetime.datetime.strftime(Stats.start, '%Y-%m-%d %H:%M:%S') +
            " - " +
            datetime.datetime.strftime(Stats.end, '%Y-%m-%d %H:%M:%S'))
    f.write("\n============SWMFU============\n")
    for M in ['mlat', 'mlt']:
        f.write(M.upper() + ':\n')
        for s in statistics:
            f.write(s + ': ' + str(round(Stats.swmfu[M].stats[s], 3)) + '\n')
    f.write("============SWMFL============\n")
    for M in ['mlat', 'mlt']:
        f.write(M.upper() + ':\n')
        for s in statistics:
            f.write(s + ': ' + str(round(Stats.swmfl[M].stats[s], 3)) + '\n')
