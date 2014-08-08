#!/usr/bin/python2

import sys
from bin.pysgpp import Grid
import GridImageFormatter

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

i = 0
for grid in grids[-5:-1]:
    formatter.serializeToFile(grid, "graph" + str(i) + ".png")
    i += 1
