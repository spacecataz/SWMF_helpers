#!/usr/bin/env python3
'''
Parse the LFS Quota command and print it in a
human-readable way.

TO MAKE THIS WORK FOR YOU: open this file and change
two variables, the "user" and the "drive", such that
they match your user name and the drive that you
are inquiring about.

DTW 2009
'''

import os

user = 'dwelling'
drive= '/nobackupp17'

# Execute lfsquota and obtain relevant data.
f=os.popen('lfs quota -u %s %s' % (user, drive))
lines=f.readlines()
for line in lines:
    if drive in line: break
f.close()

for char in ['[', ']', '*']:
    line = line.replace(char, '')

parts=line.split()
MemUsed=int(parts[1])/1000000.0
MemMax=int(parts[2])/1000000.0
nUsed=int(parts[5])
nMax=int(parts[6])

print("Used %06.2f of %06.2f Gbs (%06.2f%%) and %i of %i files (%06.2f%%)\n" %
      (MemUsed, MemMax, MemUsed/MemMax*100.0, nUsed, nMax, float(nUsed)/float(nMax)*100.0))

