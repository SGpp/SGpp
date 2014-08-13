#!/usr/bin/python2

import sys
import math
from bin.pysgpp import Grid
import GridImageFormatter
NUM_GRIDS = 10

fgrid = open(sys.argv[1], 'r')
grid_dir = sys.argv[2]

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
    print("Create Grid: " + str(i))
    formatter.serializeToFile(grids[i], grid_dir + '/' + '{:04d}'.format(i) + ".grid.png")
    k += 1
    i += interval

if i != len(grids)-1:
    i = len(grids)-1
    print("Create Grid: " + str(i))
    formatter.serializeToFile(grids[i], grid_dir + '/' + '{:04d}'.format(i) + ".grid.png")
