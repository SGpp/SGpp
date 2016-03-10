#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

# import modules
from pysgpp import DataVector, Grid, createOperationHierarchisation, \
                   createOperationQuadrature, OperationQuadratureMC

# the standard parabola (arbitrary-dimensional)
def f(x):
     res = 1.0
     for i in range(len(x)):
          res *= 4.0*x[i]*(1.0-x[i])
     return res

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:        {}".format(dim)

# create regular grid, level 3
level = 3
gridGen = grid.getGenerator()
gridGen.regular(level)
print "number of grid points: {}".format(gridStorage.getSize())

# create coefficient vector
alpha = DataVector(gridStorage.getSize())
for i in xrange(gridStorage.getSize()):
    gp = gridStorage.get(i)
    alpha[i] = f((gp.getCoord(0), gp.getCoord(1)))
createOperationHierarchisation(grid).doHierarchisation(alpha)

# direct quadrature
opQ = createOperationQuadrature(grid)
res = opQ.doQuadrature(alpha)
print "exact integral value:  {}".format(res)

# Monte Carlo quadrature using 100000 paths
opMC = OperationQuadratureMC(grid, 100000)
res = opMC.doQuadrature(alpha)
print "Monte Carlo value:     {:.6f}".format(res)
res = opMC.doQuadrature(alpha)
print "Monte Carlo value:     {:.6f}".format(res)

# Monte Carlo quadrature of a function
res = opMC.doQuadratureFunc(f)
print "MC value (f):          {:.6f}".format(res)

# Monte Carlo quadrature of error
res = opMC.doQuadratureL2Error(f, alpha)
print "MC L2-error (f-u)      {:.7f}".format(res)