# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python

# import modules
import sys
# append trunk/bin to search path for modules
sys.path.append('../lib/pysgpp')
from pysgpp import *

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:         %d" % (dim)

# create regular grid, level 3
level = 3
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "Start: number of grid points:  %d" % (gridStorage.size())

# definition of function to interpolate - nonsymmetric(!)
f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1*x1
# create coefficient vector
alpha = DataVector(gridStorage.size())

# now refine adaptively 5 times
for refnum in range(5):
    # set function values in alpha
    for i in xrange(gridStorage.size()):
        gp = gridStorage.get(i)
        alpha[i] = f(gp.abs(0), gp.abs(1))

    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)

    # refine a single grid point each time
    gridGen.refine(SurplusRefinementFunctor(alpha, 1))
    print "Refinement step %d, new grid size: %d" % (refnum+1, gridStorage.size())

    # extend alpha vector (new entries uninitialized)
    alpha.resize(gridStorage.size())