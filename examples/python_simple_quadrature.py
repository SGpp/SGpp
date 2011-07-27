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

# create coefficient vector
alpha = DataVector(gridStorage.size())
alpha.setAll(1.0)

# direct quadrature
opQ = createOperationQuadrature(grid)
res = opQ.doQuadrature(alpha)
print "exact integral value:", res

# Monte Carlo quadrature using 100000 paths
opMC = OperationQuadratureMC(grid, 100000)
res = opMC.doQuadrature(alpha)
print "Monte Carlo value:   ", res
res = opMC.doQuadrature(alpha)
print "Monte Carlo value:   ", res
