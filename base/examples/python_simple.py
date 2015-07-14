#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# import modules
from pysgpp import DataVector, Grid, createOperationHierarchisation, createOperationEval

# create a two-dimensional piecewise bilinear grid
dim = 2
grid = Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:         {}".format(gridStorage.dim())

# create regular grid, level 3
level = 3
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "number of grid points:  {}".format(gridStorage.size())

# create coefficient vector
alpha = DataVector(gridStorage.size())
alpha.setAll(0.0)
print "length of alpha vector: {}".format(len(alpha))

# set function values in alpha
f = lambda x0, x1: 16.0 * (x0-1.0)*x0 * (x1-1.0)*x1
for i in xrange(gridStorage.size()):
    gp = gridStorage.get(i)
    alpha[i] = f(gp.getCoord(0), gp.getCoord(1))
print "alpha before hierarchization: {}".format(alpha)

# hierarchize
createOperationHierarchisation(grid).doHierarchisation(alpha)
print "alpha after hierarchization:  {}".format(alpha)

# evaluate
p = DataVector(dim)
p[0] = 0.52
p[1] = 0.73
opEval = createOperationEval(grid)
print "u(0.52, 0.73) = {}".format(opEval.eval(alpha, p))
print "f(0.52, 0.73) = {}".format(f(p[0], p[1]))
