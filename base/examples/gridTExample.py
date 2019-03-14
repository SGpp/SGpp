#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_gridTExample_py Generalised Sparse Grids
## This example creates a generalised grid. The first command line argument
## is the number of dimensions, the second is the level, and the third is
## the parameter T. It then prints out the grid size.

from pysgpp import *
import sys

if len(sys.argv) != 4:
	print("The number of arguments is wrong, but the answer is 42")
else:
	dimensions = int(sys.argv[1])
	level = int(sys.argv[2])
	T = float(sys.argv[3])

	grid = Grid.createModLinearGrid(dimensions)
	generator = grid.getGenerator()
	generator.regular(level, T)
	print(grid.getSize())
