#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_quadrature_py Quadrature in Python
##
## The following example shows how to integrate in SG++, using both
## direct integration of a sparse grid function and the use of
## Monte Carlo integration.
##
## As in the \ref example_tutorial_py example, we deal with the function
## \f[
##   f\colon [0, 1]^2 \to \mathbb{R},\quad
##   f(x_0, x_1) := 16 (x_0 - 1) x_0 (x_1 - 1) x_1
## \f]
## which we first interpolate. We then integrate the interpolant, then
## the function itself using 100000 Monte Carlo points, and we then
## compute the L2-error.
##
## For instructions on how to run the example, please see \ref installation.
##
## The function, which sgpp::base::OperationQuadratureMC takes, has one parameter,
## a sequence (C++ provides a tuple) with the coordinates of the grid point \f$\in [0,1]^d\f$.
##
## This example can be found in the file quadrature.py 

# import pysgpp library
import pysgpp

# the standard parabola (arbitrary-dimensional)
def f(x):
  res = 1.0
  for i in range(len(x)):
    res *= 4.0*x[i]*(1.0-x[i])
  return res

# a pyramid-like shape (arbitrary-dimensional)
def g(x):
  res = 1.0
  for i in range(len(x)):
    res *= 2.0*min(x[i], 1.0-x[i])
  return res


## Create a two-dimensional piecewise bi-linear grid with level 3
dim = 2
grid = pysgpp.Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print("dimensionality:        {}".format(dim))

# create regular grid, level 3
level = 3
gridGen = grid.getGenerator()
gridGen.regular(level)
print("number of grid points: {}".format(gridStorage.getSize()))


## Calculate the surplus vector alpha for the interpolant of \f$
## f(x)\f$.  Since the function can be evaluated at any
## point. Hence. we simply evaluate it at the coordinates of the
## grid points to obtain the nodal values. Then we use
## hierarchization to obtain the surplus value.

# create coefficient vector
alpha = pysgpp.DataVector(gridStorage.getSize())
for i in range(gridStorage.getSize()):
  gp = gridStorage.getPoint(i)
  p = tuple([gp.getStandardCoordinate(j) for j in range(dim)])
  alpha[i] = f(p)

pysgpp.createOperationHierarchisation(grid).doHierarchisation(alpha)


## Now we compute and compare the quadrature using four different methods available in SG++.

# direct quadrature
opQ = pysgpp.createOperationQuadrature(grid)
res = opQ.doQuadrature(alpha)
print("exact integral value:  {}".format(res))

# Monte Carlo quadrature using 100000 paths
opMC = pysgpp.OperationQuadratureMC(grid, 100000)
res = opMC.doQuadrature(alpha)
print("Monte Carlo value:     {:.6f}".format(res))
res = opMC.doQuadrature(alpha)
print("Monte Carlo value:     {:.6f}".format(res))


# Monte Carlo quadrature of a standard parabola
res = opMC.doQuadratureFunc(f)
print("MC value (f):          {:.6f}".format(res))

# Monte Carlo quadrature of error
res = opMC.doQuadratureL2Error(f, alpha)
print("MC L2-error (f-u)      {:.7f}".format(res))


# Monte Carlo quadrature of a piramidal function
res = opMC.doQuadratureFunc(g)
print( "MC value (g):          {:.6f}".format(res))

# Monte Carlo quadrature of error
res = opMC.doQuadratureL2Error(g, alpha)
print( "MC L2-error (g-u)      {:.7f}".format(res))

## This results in an output similar to:
## \verbinclude quadrature.output.txt
