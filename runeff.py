#!/usr/bin/env python3
'''
Examine an SWMF logfile and determine the simulation efficiency as a percent
real time.  Basic values are written to screen; a plot in PNG format can be
saved to the pwd in order to examine the time history of the efficiency.
Efficiency is defined as the simulation time normalized by the CPU time.

Usage:

run_eff.py logfile

Options: 
-h    Print this help information.
-p    Save PNG formatted plot of efficiency and save in PWD.
-i    Interactive mode: create interactive plot.

This script requires the python module Numpy.  Though not a standard package, it
is found ubiquitously throughout the scientific community.  Plotting requires
Matplotlib; again common in scientific python environments.
'''

import os
import sys
import numpy as np

doPlot=False
doInteract=False
logfile=None

# Parse arguments.
if len(sys.argv)<=1:
    print('Correct usage: run_eff.py logfile')
    print('Use -h to print help.')
    exit()
for arg in sys.argv[1:]:
    if arg[0:2].lower()=='-h':
        print(__doc__)
        exit()
    elif arg=='-p':
        doPlot=True
    elif arg=='-i':
        doInteract=True
    else:
        logfile=arg

# Extract data points from logfile.
f=os.popen("grep '^Progress' %s"%logfile)
lines=f.readlines()
f.close()

nProg=len(lines)
cpu_t=np.zeros(nProg); run_t=np.zeros(nProg)

if nProg==0 or nProg==1:
    print("ERROR: No valid lines found.")
    print("Are you sure %s is a valid SWMF printout?" % logfile)
    exit()

for i, l in enumerate(lines):
    parts=l.replace(':', ' ').split()
    run_t[i]=float(parts[3])
    cpu_t[i]=float(parts[7])

# Shift times to include this session ONLY.
run_t=run_t-run_t[0]
cpu_t=cpu_t-cpu_t[0]
# Remove t=0 points.
good = (run_t>0)&(cpu_t>0)
run_t=run_t[good]
cpu_t=cpu_t[good]
# Efficiency; no divide-by-zero or steady-state steps.
eff=run_t/cpu_t
eff_inst = (run_t[1:] - run_t[:-1])/(cpu_t[1:] - cpu_t[:-1])
if eff.size==0:
    print("ERROR: No simulation progression detected.")
    print("Was this a steady state run?")
    exit()

# Find column number for nCPU info:
f=os.popen(f"grep 'stride' {logfile}")
lines=f.readlines()
f.close()
parts = lines[-1].split()
iCol = parts.index('nproc') - len(parts)

# Extract number of CPUs.
f=os.popen(f"grep 'CON SWMF' {logfile}")
lines=f.readlines()
f.close()
parts = lines[-1].replace('#','').split()
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
print(("Efficiency is %8.2E per CPU."   % EffCpu).center(55))
print((f"This simulation is using {nCpus:.0f} cores.").center(55))
print(("%i CPUs required for Real Time." % nCpuRT).center(55))

#print " Average Efficiency=%6.3fX Real Time" % (eff.mean())
#print " Maximum Efficiency=%6.3fX Real Time" % (eff.max())
#print " Minimum Efficiency=%6.3fX Real Time" % (eff.min())


if doPlot:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from scipy.interpolate import interp1d

    # Create custom ticks that have sim and cpu time.
    cpu2sim = interp1d(cpu_t/3600.0, run_t/3600.0, bounds_error=False, 
                       fill_value=0)
    def cust_tick(x, pos):
        sim_time = cpu2sim(x)
        return '%5.2f\n%5.2f' % (x, sim_time)

    f=plt.figure()
    f.subplots_adjust(top=0.9, bottom=0.15, left=0.1, right=0.95)
    a1=f.add_subplot(111)
    a1.plot(cpu_t/3600.0, eff, 'b.-', label='Cumulative Eff.', zorder=100)
    #a1.plot(cpu_t[1:]/3600.0, eff_inst, 'g-', label='_nolegend_',
    #        zorder=10, alpha=0.25)
    a1.plot(cpu_t[1:]/3600.0, eff_inst, 'g.', label='Instantaneous Eff.',
            zorder=10, ms=1.2)
    a1.hlines(eff[-1], cpu_t[run_t>0][0]/3600., cpu_t[run_t>0][-1]/3600., 
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
    
    if doInteract:
        plt.show()
    else:
        f.savefig('efficiency.png')
