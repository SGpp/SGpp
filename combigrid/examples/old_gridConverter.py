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

import numpy as np
import matplotlib.pyplot as plt

from pysgpp import Grid, fullStorageLevels, allStorageLevels, toHashGridStorage
from pysgpp.pysgpp_swig import CombigridOperation, multiFunc, DataVector, \
    CombigridMultiOperation, DataMatrix, SurplusRefinementFunctor
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSG1d
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hierarchize, \
    evalSGFunction

numDims = 1
level = 2
refnums = 0
n = 100

## Then, we define the problem: we want to interpolate the function f
def f(x):
    return np.prod([4 * xi * (1 - xi) for xi in x.array()])

grid = Grid.createLinearGrid(numDims)
gridGen = grid.getGenerator()
gridGen.regular(level)
regularGridSize = grid.getSize()

gs = grid.getStorage()

x = np.linspace(0, 1, n)
x.resize(1, n)

# x = np.random.rand(numDims, n)
parameters = DataMatrix(x)

nodalValues = np.ndarray(gs.getSize())
p = DataVector(numDims)
for i in xrange(gs.getSize()):
    gs.getCoordinates(gs.getPoint(i), p)
    nodalValues[i] = f(p)

alpha = hierarchize(grid, nodalValues)

for _ in xrange(refnums):
    # refine a single grid point each time
    alpha_vec = DataVector(alpha)
    gridGen.refine(SurplusRefinementFunctor(alpha_vec, 1))

    # extend alpha vector (new entries uninitialized)
    nodalValues.resize(gs.getSize())

    # set function values in alpha
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        nodalValues[i] = f(p)

    # hierarchize
    alpha = hierarchize(grid, nodalValues)

y_sg = evalSGFunction(grid, alpha, x.T)

treeStorage = allStorageLevels(grid.getStorage())
# treeStorage = fullStorageLevels(grid.getStorage())
opt = CombigridMultiOperation.createExpUniformLinearInterpolation(numDims, multiFunc(f))
opt.setParameters(parameters)
opt.getLevelManager().addLevelsFromStructure(treeStorage)
y = opt.getResult().array()

# print opt.getLevelManager().getSerializedLevelStructure()
# if refnums == 0:
#     assert np.sum((y - y_sg) ** 2) < 1e-14

newGrid = Grid.createLinearGrid(numDims)
treeStorage = opt.getLevelManager().getLevelStructure()
toHashGridStorage(treeStorage, newGrid.getStorage())

newGrid.getStorage().recalcLeafProperty()
newGs = newGrid.getStorage()
nodalValues.resize(newGs.getSize())
for i in xrange(newGs.getSize()):
    newGs.getCoordinates(newGs.getPoint(i), p)
    nodalValues[i] = f(p)
newAlpha = hierarchize(newGrid, nodalValues)

x = x.flatten()
ixs = np.argsort(x)
plt.figure()
plt.plot(x[ixs], y_sg[ixs], label="regular")
plt.plot(x[ixs], y[ixs], label="full")
plt.legend()
plt.show()


# if refnums == 0:
#     y_sg_new = evalSGFunction(newGrid, newAlpha, x.T)
#     assert np.sum((y - y_sg_new) ** 2) < 1e-14
# else:
#     assert newGrid.getSize() == regularGridSize
