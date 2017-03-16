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
from argparse import ArgumentParser
import matplotlib.pyplot as plt

from pysgpp import CombigridOperation, multiFunc, DataVector, \
    CombigridMultiOperation, DataMatrix, SurplusRefinementFunctor, \
    Grid, fullStorageLevels, allStorageLevels, toHashGridStorage, \
    createOperationHierarchisation, createOperationMultipleEval

## Then, we define the problem: we want to interpolate the function f
def f(x):
    return np.prod([4 * xi * (1 - xi) for xi in x.array()])


def interpolate(grid, f):
    ## We interpolate now the function on the sparse grid in the
    ## hierarchical version and...
    gs = grid.getStorage()
    alpha = DataVector(gs.getSize())
    p = DataVector(gs.getDimension())
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        alpha[i] = f(p)
    createOperationHierarchisation(grid).doHierarchisation(alpha)
    return alpha


def refineGrid(grid, alpha, f, refnums):
    gs = grid.getStorage()
    gridGen = grid.getGenerator()
    x = DataVector(gs.getDimension())
    for _ in xrange(refnums):
        # refine a single grid point each time
        gridGen.refine(SurplusRefinementFunctor(alpha, 1))
    
        # extend alpha vector (new entries uninitialized)
        alpha.resize(gs.getSize())
    
        # set function values in alpha
        for i in xrange(gs.getSize()):
            gs.getCoordinates(gs.getPoint(i), x)
            alpha[i] = f(x)
    
        # hierarchize
        createOperationHierarchisation(grid).doHierarchisation(alpha)


def regularGridToRegularGrid(numDims, level, f, n=1000):
    # # We generate a iid of uniform samples, which we are going to use to
    # # validate the grid conversion
    x = np.random.rand(n, numDims)
    parameters = DataMatrix(x)
    
    ## We create a regular sparse grid as usual...
    grid = Grid.createLinearGrid(numDims)
    grid.getGenerator().regular(level)
    alpha = interpolate(grid, f)
    
    # # We apply now both methods of the grid conversion
    treeStorage_all = allStorageLevels(grid.getStorage())
    treeStorage_full = fullStorageLevels(grid.getStorage())
    
    # # Note, that we do the conversion just based on the grid points. With
    # # this approach you can easily use different basis function types on
    # # the same grid. We initialize the CombigridOperation on a grid that
    # # spans the same function space as the original hierarchical sparse
    # # grid: hat basis on an equidistant grids without boundary points.
    func = multiFunc(f)
    opt_full = CombigridMultiOperation.createExpUniformLinearInterpolation(numDims, func)
    opt_all = CombigridMultiOperation.createExpUniformLinearInterpolation(numDims, func)
    
    # # The CombigridOperation expects the points at which you want to
    # # evaluate the interpolant as DataMatrix with the shape (numDims x
    # # numSamples). We initialize the levels for both operations.
    
    parameters.transpose()
    opt_full.setParameters(parameters)
    opt_full.getLevelManager().addLevelsFromStructure(treeStorage_full)
    opt_all.setParameters(parameters)
    opt_all.getLevelManager().addLevelsFromStructure(treeStorage_all)
    parameters.transpose()
    
    import idpb; ipdb.set_trace()
    print "-" * 80
    print "just full levels:"
    print opt_full.getLevelManager().getSerializedLevelStructure()
    print "-" * 80
    print "all levels:"
    print opt_all.getLevelManager().getSerializedLevelStructure()
    print "-" * 80

    grid_full = Grid.createLinearGrid(numDims)
    gs_full = grid_full.getStorage()
    treeStorage_full = opt_full.getLevelManager().getLevelStructure()
    toHashGridStorage(treeStorage_full, gs_full)

    grid_all = Grid.createLinearGrid(numDims)
    gs_all = grid_all.getStorage()
    treeStorage_all = opt_all.getLevelManager().getLevelStructure()
    toHashGridStorage(treeStorage_all, gs_all)

    alpha_full = interpolate(grid_full, f)
    alpha_all = interpolate(grid_all, f)

    # # evaluate the surrogate functions
    y_sg_regular = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid, parameters).eval(alpha, y_sg_regular)
    
    y_sg_all = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid_all, parameters).eval(alpha_all, y_sg_all)
    
    y_sg_full = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid_full, parameters).eval(alpha_full, y_sg_full)

    y_ct_all = opt_all.getResult()
    y_ct_full = opt_full.getResult()

    # # tests
    y_sg_regular = y_sg_regular.array().flatten()
    y_ct_all = y_ct_all.array().flatten()
    y_ct_full = y_ct_full.array().flatten()
    y_sg_all = y_sg_all.array().flatten()
    y_sg_full = y_sg_full.array().flatten()

    if numDims == 1:
        x = x.flatten()
        ixs = np.argsort(x)
        plt.figure()
        plt.plot(x[ixs], y_sg_regular[ixs], label="regular")
        plt.plot(x[ixs], y_sg_all[ixs], label="sg all")
        plt.plot(x[ixs], y_ct_full[ixs], label="full")
        plt.plot(x[ixs], y_ct_all[ixs], label="all")
        plt.legend()
        plt.show()

#     assert np.sum((y_ct_full - y_ct_all) ** 2) < 1e-14
#     assert np.sum((y_ct_full - y_sg_regular) ** 2) < 1e-14
    assert np.sum((y_sg_regular - y_sg_all) ** 2) < 1e-14
    assert np.sum((y_sg_regular - y_sg_full) ** 2) < 1e-14
    assert grid_full.getSize() == grid.getSize()
    assert grid_all.getSize() == grid.getSize()


def adaptiveGridToRegularGrid(numDims, level, refnums, f, n=1000):
    # # We generate a iid of uniform samples, which we are going to use to
    # # validate the grid conversion
    x = np.random.rand(n, numDims)
    parameters = DataMatrix(x)

    # # We create a regular sparse grid as usual...
    grid = Grid.createLinearGrid(numDims)
    grid.getGenerator().regular(level)
    alpha = interpolate(grid, f)

    grid_adaptive = grid.clone()
    alpha_adaptive = DataVector(alpha)
    refineGrid(grid_adaptive, alpha_adaptive, f, refnums)

    # # We apply now both methods of the grid conversion
    treeStorage_all = allStorageLevels(grid_adaptive.getStorage())
    treeStorage_full = fullStorageLevels(grid_adaptive.getStorage())

    # # Note, that we do the conversion just based on the grid points. With
    # # this approach you can easily use different basis function types on
    # # the same grid. We initialize the CombigridOperation on a grid that
    # # spans the same function space as the original hierarchical sparse
    # # grid: hat basis on an equidistant grids without boundary points.
    func = multiFunc(f)
    opt_full = CombigridMultiOperation.createExpUniformLinearInterpolation(numDims, func)
    opt_all = CombigridMultiOperation.createExpUniformLinearInterpolation(numDims, func)

    # # The CombigridOperation expects the points at which you want to
    # # evaluate the interpolant as DataMatrix with the shape (numDims x
    # # numSamples). We initialize the levels for both operations.

    parameters.transpose()
    opt_full.setParameters(parameters)
    opt_full.getLevelManager().addLevelsFromStructure(treeStorage_full)
    opt_all.setParameters(parameters)
    opt_all.getLevelManager().addLevelsFromStructure(treeStorage_all)
    parameters.transpose()

    print "-" * 80
    print "just full levels:"
    print opt_full.getLevelManager().getSerializedLevelStructure()
    print "-" * 80
    print "all levels:"
    print opt_all.getLevelManager().getSerializedLevelStructure()
    print "-" * 80

    grid_full = Grid.createLinearGrid(numDims)
    treeStorage_full = opt_full.getLevelManager().getLevelStructure()
    toHashGridStorage(treeStorage_full, grid_full.getStorage())

    grid_all = Grid.createLinearGrid(numDims)
    treeStorage_all = opt_all.getLevelManager().getLevelStructure()
    toHashGridStorage(treeStorage_all, grid_all.getStorage())

    alpha_full = interpolate(grid_full, f)
    alpha_all = interpolate(grid_all, f)

    # # evaluate the surrogate functions
    y_sg_regular = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid, parameters).eval(alpha, y_sg_regular)

    y_sg_adaptive = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid_adaptive, parameters).eval(alpha_adaptive, y_sg_adaptive)

    y_sg_all = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid_all, parameters).eval(alpha_all, y_sg_all)

    y_sg_full = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid_full, parameters).eval(alpha_full, y_sg_full)

    y_ct_all = opt_all.getResult()
    y_ct_full = opt_full.getResult()

    # # tests
    y_sg_regular = y_sg_regular.array().flatten()
    y_sg_adaptive = y_sg_adaptive.array().flatten()
    y_ct_all = y_ct_all.array().flatten()
    y_ct_full = y_ct_full.array().flatten()
    y_sg_all = y_sg_all.array().flatten()
    y_sg_full = y_sg_full.array().flatten()

    if numDims == 1:
        x = x.flatten()
        ixs = np.argsort(x)
        plt.figure()
#         plt.plot(x[ixs], y_sg_regular[ixs], label="regular")
        plt.plot(x[ixs], y_sg_adaptive[ixs], label="adaptive")
#         plt.plot(x[ixs], y_sg_full[ixs], label="sg full")
#         plt.plot(x[ixs], y_sg_all[ixs], label="sg all")
        plt.plot(x[ixs], y_ct_full[ixs], label="full")
#         plt.plot(x[ixs], y_ct_all[ixs], label="all")
        plt.legend()
        plt.show()

#     assert np.sum((y_ct_full - y_ct_all) ** 2) < 1e-14
#     assert np.sum((y_ct_full - y_sg_regular) ** 2) < 1e-14
#     assert np.sum((y_sg_regular - y_sg_all) ** 2) > 1e-14
#     assert np.sum((y_sg_regular - y_sg_full) ** 2) < 1e-14

    assert grid_adaptive.getSize() > grid.getSize()
    assert grid_full.getSize() <= grid_adaptive.getSize()
    assert grid_all.getSize() >= grid.getSize()

    if grid_full.getSize() == grid.getSize():
        assert np.sum((y_ct_full - y_sg_regular) ** 2) < 1e-14


if __name__ == '__main__':
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--numDims', default=4, type=int, help='number of dimensions')
    parser.add_argument('--level', default=4, type=int, help='sparse grid level')
    parser.add_argument('--refnums', default=0, type=int, help='number of refinement steps')
    args = parser.parse_args()

    if args.refnums == 0:
        regularGridToRegularGrid(args.numDims,
                                 args.level,
                                 f)
    else:
        adaptiveGridToRegularGrid(args.numDims,
                                  args.level,
                                  args.refnums,
                                  f)
