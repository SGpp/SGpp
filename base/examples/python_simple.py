# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python

# import modules
import sys
# append trunk/bin to search path for modules
sys.path.append('../lib/pysgpp')
from pysgpp.base import DataVector, Grid, createOperationHierarchisation, createOperationEval

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:         %d" % (gridStorage.dim())

# create regular grid, level 3
level = 3
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "number of grid points:  %d" % (gridStorage.size())

# create coefficient vector
alpha = DataVector(gridStorage.size())
alpha.setAll(0.0)
print "length of alpha-vector: %d" % (len(alpha))

# set function values in alpha
f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1
for i in xrange(gridStorage.size()):
    gp = gridStorage.get(i)
    alpha[i] = f(gp.abs(0), gp.abs(1))
print alpha

# hierarchize
createOperationHierarchisation(grid).doHierarchisation(alpha)
print alpha

# evaluate
p = DataVector(dim)
p[0] = 0.52
p[1] = 0.73
opEval = createOperationEval(grid)
print "u(0.52, 0.73) =", opEval.eval(alpha, p)