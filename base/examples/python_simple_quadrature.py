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

# the standard parabola (arbitrary-dimensional)
def f(x):
     res = 1.0
     for i in range(len(x)):
          res *= 4*x[i]*(1-x[i])
     return res

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:        %d" % (dim)

# create regular grid, level 3
level = 3
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "number of grid points: %d" % (gridStorage.size())

# create coefficient vector
alpha = DataVector(gridStorage.size())
for i in xrange(gridStorage.size()):
    gp = gridStorage.get(i)
    alpha[i] = f((gp.getCoord(0), gp.getCoord(1)))
createOperationHierarchisation(grid).doHierarchisation(alpha)

# direct quadrature
opQ = createOperationQuadrature(grid)
res = opQ.doQuadrature(alpha)
print "exact integral value: ", res

# Monte Carlo quadrature using 100000 paths
opMC = OperationQuadratureMC(grid, 100000)
res = opMC.doQuadrature(alpha)
print "Monte Carlo value:    ", res
res = opMC.doQuadrature(alpha)
print "Monte Carlo value:    ", res

# Monte Carlo quadrature of a function
res = opMC.doQuadratureFunc(f)
print "MC value (f):         ", res

# Monte Carlo quadrature of error
res = opMC.doQuadratureL2Error(f, alpha)
print "MC L2-error (f-u)     ", res