#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_refinement_py refinement.py
##
## The following example interpolates the (non-symmetric) function
## \f[
##   f\colon [0, 1]^2 \to \mathbb{R},\quad
##   f(x_0, x_1) := 16 (x_0 - 1) x_0 (x_1 - 1) x_1
## \f]
## - starting from a regular grid with level 2 and
## - refining 5 times one grid point each.
##
## The number of grid points is printed in each iteration.
## After refinement, the surplusses have to be set for all new grid
## points, i.e., the alpha-Vector has to be extended.
##
## For instructions on how to run the example, please see \ref installation.

# import pysgpp library
import pysgpp

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = pysgpp.Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:                   {}".format(dim)

# create regular grid, level 3
level = 3
gridGen = grid.getGenerator()
gridGen.regular(level)
print "number of initial grid points:    {}".format(gridStorage.getSize())

# definition of function to interpolate - nonsymmetric(!)
f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1*x1
# create coefficient vector
alpha = pysgpp.DataVector(gridStorage.getSize())
print "length of alpha vector:           {}".format(alpha.getSize())


# now refine adaptively 5 times
for refnum in range(5):
    # set function values in alpha
    for i in xrange(gridStorage.getSize()):
        gp = gridStorage.getPoint(i)
        alpha[i] = f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1))

    # hierarchize
    pysgpp.createOperationHierarchisation(grid).doHierarchisation(alpha)

    # refine a single grid point each time
    gridGen.refine(pysgpp.SurplusRefinementFunctor(alpha, 1))
    print "refinement step {}, new grid size: {}".format(refnum+1, gridStorage.getSize())

    # extend alpha vector (new entries uninitialized)
    alpha.resize(gridStorage.getSize())

## This results in the following output:
## \verbinclude refinement.output.txt
##
## There are clearly more efficient approaches than to set the function
## values for all grid points and to hierarchize the whole grid each
## time. But this works even where no efficient alternatives are
## available and suffices for demonstration purposes.
##
## This use of the SurplusRefinementFunctor takes as arguments the
## coefficient vector (it doesn't have to be the coefficient vector, it
## could be something modified!) and the number of grid points to refine
## (if available). It bases its refinement decision on the absolute
## values of the vector's entries, choosing the largest ones. Other
## refinement functors are available or can be implemented.
