#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_refinement_py refinement.py
##
## Here we demonstrate how to refine a grid. As a refinement indicator, we take the surpluses
## of the grid points directly.
## We start with a regular sparse grid of level 3 with linear basis functions and refine five
## times. In each refinement step, we refine the grid point with the highest absolute surplus.
##
## The following example interpolates the (non-symmetric) function
## \f[
##   f\colon [0, 1]^2 \to \mathbb{R},\quad
##   f(x_0, x_1) := 16 (x_0 - 1) x_0 (x_1 - 1) x_1
## \f]
##
## The number of grid points is printed in each iteration.
## After refinement, the surplusses have to be set for all new grid
## points, i.e., the alpha-Vector has to be extended.

# import pysgpp library
import pysgpp

## create a two-dimensional piecewise bi-linear grid
dim = 2
grid = pysgpp.Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print("dimensionality:                   {}".format(dim))

## create regular sparse grid, level 3
level = 3
gridGen = grid.getGenerator()
gridGen.regular(level)
print("number of initial grid points:    {}".format(gridStorage.getSize()))

## function to interpolate. This is a two-dimensional parabola. - nonsymmetric(!)
f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1*x1
## create coefficient vector with size corresponding to the grid size.
## Initially, all the values are set to zero.
alpha = pysgpp.DataVector(gridStorage.getSize())
print("length of alpha vector:           {}".format(alpha.getSize()))


# Obtain function values and refine adaptively 5 times
for refnum in range(5):
    # set function values in alpha
    for i in range(gridStorage.getSize()):
        gp = gridStorage.getPoint(i)
        alpha[i] = f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1))

    ## Each time, we have to hierarchize the grid again, because in the previous interation,
    ## new grid points have been added.
    pysgpp.createOperationHierarchisation(grid).doHierarchisation(alpha)

    ## Refine a single grid point each time.
    ## The SurplusRefinementFunctor chooses the grid point with the highest absolute surplus.
    ## Refining the point means, that all children of this point (if not already present) are
    ## added to the grid. Also all missing parents are added (recursively).
    gridGen.refine(pysgpp.SurplusRefinementFunctor(alpha, 1))
    print("refinement step {}, new grid size: {}".format(refnum+1, gridStorage.getSize()))

    ## Extend alpha vector (new entries uninitialized). Note that right now, the surplus vector
    ## has the correct size again, but the values of the new points are set to zero. The correct
    ## surplus values will be inserted in the next iteration of the refinement loop.
    alpha.resizeZero(gridStorage.getSize())

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
