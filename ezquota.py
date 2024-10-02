#!/usr/bin/env python3
'''
Parse the LFS Quota command and print it in a
human-readable way.

DTW 2009
'''

import os

# First, grab current user:
user = os.getlogin()
drive = '/nobackup'

# Execute lfsquota and obtain relevant data.
f = os.popen('lfs quota -u {} {}'.format(user, drive))
lines = f.readlines()
for line in lines:
    if drive in line:
        break
f.close()

for char in ['[', ']', '*']:
    line = line.replace(char, '')

parts = line.split()
MemUsed = int(parts[1])/1000000.0
MemMax = int(parts[2])/1000000.0
nUsed = int(parts[5])
nMax = int(parts[6])

print(f"Used {MemUsed:06.2f} of {MemMax:06.2f} Gbs " +
      f"({MemUsed/MemMax*100.0:06.2f}%) and {nUsed:d} of {nMax:d} " +
      f"files ({float(nUsed)/float(nMax)*100.0:06.2f}%)\n")

