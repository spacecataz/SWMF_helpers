#!/usr/bin/env python3
'''
Examine an SWMF logfile and determine the simulation efficiency as a percent
real time.  Basic values are written to screen; a plot in PNG format can be
saved to the pwd in order to examine the time history of the efficiency.
Efficiency is defined as the simulation time normalized by the CPU time.

This script requires the python module Numpy.  Though not a standard package,
it is found ubiquitously throughout the scientific community.  Plotting
requires Matplotlib; again common in scientific python environments.
'''

import os
import numpy as np
from argparse import ArgumentParser

# Set and parse arguments:
parser = ArgumentParser(description=__doc__)
parser.add_argument("-p", "--plot",  action='store_true',
                    help="Save PNG formatted plot of efficiency and save.")
parser.add_argument("-i", "--interactive", action="store_true",
                    help="Interactive mode: create interactive plot.")
parser.add_argument("logfile", nargs=1, help="Run log file to scan.")
parser.add_argument("-d", "--debug", action="store_true",
                    help="Set debug output.")
args = parser.parse_args()

# Set 1 and only 1 log file:
args.logfile = args.logfile[0]

# Extract data points from logfile.
if args.debug:
    print(f"Reading log file {args.logfile}")
with os.popen(f"grep -a '^Progress' {args.logfile}") as f:
    lines = f.readlines()

nProg = len(lines)
cpu_t = np.zeros(nProg)
run_t = np.zeros(nProg)
dates = []

if nProg == 0 or nProg == 1:
    print("ERROR: No valid lines found.")
    print(f"Are you sure {args.logfile} is a valid SWMF printout?")
    exit()

for i, l in enumerate(lines):
    parts = l.replace(':', ' ').split()
    run_t[i] = float(parts[3])
    cpu_t[i] = float(parts[7])
    dates.append(parts[-1])

# Shift times to include this session ONLY.
run_t = run_t-run_t[0]
cpu_t = cpu_t-cpu_t[0]

# Remove t=0 points.
good = (run_t > 0) & (cpu_t > 0)
run_t = run_t[good]
cpu_t = cpu_t[good]

# Efficiency; no divide-by-zero or steady-state steps.
eff = run_t/cpu_t
eff_inst = (run_t[1:] - run_t[:-1])/(cpu_t[1:] - cpu_t[:-1])
if eff.size == 0:
    print("ERROR: No simulation progression detected.")
    print("Was this a steady state run?")
    exit()

# Find column number for nCPU info:
with os.popen(f"grep -a 'stride' {args.logfile}") as f:
    lines = f.readlines()
parts = lines[-1].split()
iCol = parts.index('nproc') - len(parts)

# Extract number of CPUs.
with os.popen(f"grep -a 'CON SWMF' {args.logfile}") as f:
    lines = f.readlines()
parts = lines[-1].split()
nCpus = float(parts[iCol])

# Get "previous hour" values.
PrevLoc = cpu_t >= (cpu_t[-1]-3600.0)
MedInst = np.median(eff_inst[PrevLoc[1:]])
EffRate = 3600.0*(eff[PrevLoc][-1]-eff[PrevLoc][0]) / \
    (cpu_t[PrevLoc][-1]-cpu_t[PrevLoc][0])
EffCpu = eff[-1]/nCpus
nCpuRT = EffCpu**-1


# Write report to screen.
print("----------=========Efficiency Report=========----------")
print(("Simulated %06.2f Hrs in %06.2f Hrs (%06.3fX Real Time)" %
       (run_t[-1]/3600.0, cpu_t[-1]/3600.0, eff[-1])).center(55))
print(("Median Instantaneous Eff. in past hour = %06.3fX" %
       (np.median(eff_inst[PrevLoc[1:]]))).center(55))
if EffRate < 0:
    print(("Efficiency is DECREASING by %8.2E per hour." %
           (-1.*EffRate)).center(55))
else:
    print(("Efficiency is INCREASING by %8.2E per hour." % EffRate).center(55))
print(("Efficiency is %8.2E per CPU." % EffCpu).center(55))
print((f"This simulation is using {nCpus:.0f} cores.").center(55))
print(("%i CPUs required for Real Time." % nCpuRT).center(55))
print((f"Start datetime is {dates[0]}").center(55))
print((f"Curr. datetime is {dates[-1]}").center(55))

#  print " Average Efficiency=%6.3fX Real Time" % (eff.mean())
#  print " Maximum Efficiency=%6.3fX Real Time" % (eff.max())
#  print " Minimum Efficiency=%6.3fX Real Time" % (eff.min())


if args.plot:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from scipy.interpolate import interp1d

    # Create custom ticks that have sim and cpu time.
    cpu2sim = interp1d(cpu_t/3600.0, run_t/3600.0, bounds_error=False,
                       fill_value=0)

    def cust_tick(x, pos):
        sim_time = cpu2sim(x)
        return '%5.2f\n%5.2f' % (x, sim_time)

    f = plt.figure()
    f.subplots_adjust(top=0.9, bottom=0.15, left=0.1, right=0.95)
    a1 = f.add_subplot(111)
    a1.plot(cpu_t/3600.0, eff, 'b.-', label='Cumulative Eff.', zorder=100)
    #  a1.plot(cpu_t[1:]/3600.0, eff_inst, 'g-', label='_nolegend_',
    #        zorder=10, alpha=0.25)
    a1.plot(cpu_t[1:]/3600.0, eff_inst, 'g.', label='Instantaneous Eff.',
            zorder=10, ms=1.2)
    a1.hlines(eff[-1], cpu_t[run_t > 0][0]/3600., cpu_t[run_t > 0][-1]/3600.,
              colors='r', linestyles='dashed', label='Final Cum. Eff.',
              zorder=101, lw=2.0)
    a1.grid()
    a1.xaxis.set_major_formatter(FuncFormatter(cust_tick))
    a1.legend(loc='best')
    a1.set_xlabel('CPU Time (Hours)\nSim Tim (Hours)')
    a1.set_ylabel('Run Speed (Sim Time/CPU Time)')
    a1.set_title("Simulated %06.2f Hrs in %06.2f Hrs (%06.3fX Real Time)" %
                 (run_t[-1]/3600.0, cpu_t[-1]/3600.0, eff[-1]) +
                 "\n%i CPUs required for Real Time." % nCpuRT)

    if args.interactive:
        plt.show()
    else:
        f.savefig('efficiency.png')
