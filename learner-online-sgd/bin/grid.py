#!/usr/bin/python2

import sys
import math
from bin.pysgpp import Grid
import GridImageFormatter
NUM_GRIDS = 10


fgrid = open(sys.argv[1], 'r')
grids = []
accum = ""
for line in fgrid:
    if line != "\n":
        accum += line
    else:
        grids.append(Grid.unserialize(accum))
        accum=""

if accum != "":
    grids.append(Grid.unserialize(accum))

formatter = GridImageFormatter.GridImageFormatter()

interval = int((len(grids) - NUM_GRIDS) / (NUM_GRIDS-1)) + 1

k = 1
i = 0
while i < len(grids):
    formatter.serializeToFile(grids[i], '{:04d}'.format(k) + ".grid.png")
    k += 1
    i += interval
