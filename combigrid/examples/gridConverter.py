#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page combigrid_gridConverter_py gridConverter.py This tutorial
## contains examples on how to convert sparse grids with a
## hierarchical basis to a sparse grid defined on the combination of
## anisotropic full grids (combination technique). It just includes
## grids without boundary points.

## We distinguish between methods that convert (1) the anisotropic full
## grids to the hierarchical grids and (2) vice vera:
## (1) Converting the levels from the combination technique to the
##     hierarchical version is always possible -> toHashGridStorage
## (2) For spatially adaptive sparse grids it is possible that there
##     exist just partially filled levels. Partially filled
##     levels are not allowed in the combination technique. We,
##     therefore, distinguish two cases where we
##     (a) Add all levels where there exists at least one grid point
##         in the hierarchical version -> allStorageLevels
##     (b) Add just those levels where exist all the grid points in
##         the hierarchical version -> fullStorageLevels

## First, we import a the methods/classes we need for this example
import numpy as np
from pysgpp import CombigridOperation, multiFunc, DataVector, \
    CombigridMultiOperation, DataMatrix, SurplusRefinementFunctor, \
    Grid, fullStorageLevels, allStorageLevels, toHashGridStorage, \
    createOperationHierarchisation, createOperationMultipleEval

## Then, we define the problem: we want to interpolate the function f
def f(x):
    return np.prod([4 * xi * (1 - xi) for xi in x])

## defined on $[0, 1]^4$ with an adaptively refined (refnum-times)
## sparse grid with a hierarchical basis
numDims = 1
level = 1
refnums = 0

# # We generate a iid of uniform samples, which we are going to use to
# # validate the grid conversion
n = 1000
x = np.random.rand(numDims, n)
parameters = DataMatrix(x)

## We create a regular sparse grid as usual...
grid = Grid.createLinearGrid(numDims)
grid.getGenerator().regular(level)
gs = grid.getStorage()

## We interpolate now the function on the sparse grid in the
## hierarchical version and...
alpha = DataVector(gs.getSize())
p = DataVector(numDims)
for i in xrange(gs.getSize()):
    gs.getCoordinates(gs.getPoint(i), p)
    alpha[i] = f(p.array())
createOperationHierarchisation(grid).doHierarchisation(alpha)

# # ... we furthermore apply spatially adaptive refinement such that we
# # create additional levels, which are, however, partially empty.
grid_adaptive = grid.clone()
alpha_adaptive = DataVector(alpha)
gs_adaptive = grid_adaptive.getStorage()
gridGen_adaptive = grid_adaptive.getGenerator()
for _ in xrange(refnums):
    # refine a single grid point each time
    gridGen_adaptive.refine(SurplusRefinementFunctor(alpha_adaptive, 1))

    # extend alpha vector (new entries uninitialized)
    alpha_adaptive.resize(gs_adaptive.getSize())

    # set function values in alpha
    for i in xrange(gs_adaptive.getSize()):
        gs_adaptive.getCoordinates(gs_adaptive.getPoint(i), p)
        alpha_adaptive[i] = f(p.array())

    # hierarchize
    createOperationHierarchisation(grid_adaptive).doHierarchisation(alpha_adaptive)

# # We apply now both methods of the grid conversion
treeStorage_all = allStorageLevels(grid.getStorage())
treeStorage_full = fullStorageLevels(grid_adaptive.getStorage())

print treeStorage_all
# it = treeStorage_all.getStoredDataIterator()
# while it.isValid():
#     print it.getMultiIndex()
#     it.moveToNext()

# # Note, that we do the conversion just based on the grid points. With
# # this approach you can easily use different basis function types on
# # the same grid. We initialize the CombigridOperation on a grid that
# # spans the same function space as the original hierarchical sparse
# # grid: hat basis on an equidistant grids without boundary points.
func = multiFunc(f)
opt_full = CombigridMultiOperation.createExpUniformNoBoundaryLinearInterpolation(numDims, func)
opt_all = CombigridMultiOperation.createExpUniformNoBoundaryLinearInterpolation(numDims, func)

# # The CombigridOperation expects the points at which you want to
# # evaluate the interpolant as DataMatrix with the shape (numDims x
# # numSamples). We initialize the levels for both operations.

parameters.transpose()
opt_full.setParameters(parameters)
# opt_full.getLevelManager().addLevelsFromStructure(treeStorage_full)
opt_all.setParameters(parameters)
# opt_all.getLevelManager().addLevelsFromStructure(treeStorage_all)
parameters.transpose()

# #
# grid_new = Grid.createLinearGrid(numDims)
# gs_new = newGrid.getStorage()
# treeStorage_full = opt_full.getLevelManager().getLevelStructure()
# toHashGridStorage(treeStorage_full, gs_new)
#
# alpha_new = DataVector(gs_new.getSize())
# for i in xrange(gs_new.getSize()):
#     gs_new.getCoordinates(gs_new.getPoint(i), p)
#     alpha_new[i] = f(p.array())
# createOperationHierarchisation(newGrid).doHierarchisation(alpha_new)
#
#
# y_sg = DataVector(parameters.getNrows())
# createOperationMultipleEval(grid, parameters).eval(alpha, y_sg)
# y_ct_all = opt_all.getResult()
# y_ct_full = opt_full.getResult()
#
# y_sg_new = DataVector(parameters.getNrows())
# createOperationMultipleEval(newGrid, parameters).eval(alpha_new, y_sg_new)
#
# y_sg = y_sg.array()
# y_ct_all = y_ct_all.array()
# y_ct_full = y_ct_all.array()
# y_sg_new = y_sg_new.array()
#
# if refnums == 0:
#     assert np.sum((y_ct - y_sg) ** 2) < 1e-14
#
# assert np.sum((y_ct_full - y_sg_new) ** 2) < 1e-14
# assert newGrid.getSize() == regularGridSize
