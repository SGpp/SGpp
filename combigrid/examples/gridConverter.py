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
##     hierarchical version is always possible
##      -> convertCombigridToHierarchicalSparseGrid
## (2) For spatially adaptive sparse grids it is possible that there
##     exist just partially filled levels. Partially filled
##     levels are not allowed in the combination technique. We,
##     therefore, distinguish two cases where we
##     (a) Add all levels where there exists at least one grid point
##         in the hierarchical version
##         -> convertHierarchicalSparseGridToCombigrid with
##            conversion type GridConversionTypes_ALLSUBSPACES
##     (b) Add just those levels where exist all the grid points in
##         the hierarchical version
##         -> convertHierarchicalSparseGridToCombigrid with
##            conversion type GridConversionTypes_COMPLETESUBSPACES

## First, we import a the methods/classes we need for this example...
import numpy as np
from argparse import ArgumentParser
import matplotlib.pyplot as plt

from pysgpp import CombigridOperation, multiFunc, DataVector, \
    CombigridMultiOperation, DataMatrix, SurplusRefinementFunctor, \
    Grid, convertCombigridToHierarchicalSparseGrid, convertHierarchicalSparseGridToCombigrid, \
    GridConversionTypes_ALLSUBSPACES, GridConversionTypes_COMPLETESUBSPACES, \
    createOperationHierarchisation, createOperationMultipleEval

## ... and define we define the function we want to interpolate. It is
## a parbola, which is zero for any $x_i=0$ and $x_i=1$ and evaluates
## to $1$ for $\vec{x} = (0.5, \dots, 0.5)^T
def f(x):
    return np.prod([4 * xi * (1 - xi) for xi in x.array()])


## @section py_combigrid_grid_converter_1 Helper functions: We first define a few functions that remove boiler plate code from the actual examples.
def interpolate(grid, f):
    """
    This helper functions cmoputes the coefficients of a sparse grid
    function for a given function

    Arguments:
    grid -- Grid sparse grid from pysgpp
    f -- function to be interpolated

    Return DataVector coefficients of the sparse grid function
    """
    gs = grid.getStorage()
    alpha = DataVector(gs.getSize())
    p = DataVector(gs.getDimension())
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        alpha[i] = f(p)
    createOperationHierarchisation(grid).doHierarchisation(alpha)
    return alpha


def refineGrid(grid, alpha, f, refnums):
    """
    This function refines a sparse grid function refnum times.

    Arguments:
    grid -- Grid sparse grid from pysgpp
    alpha -- DataVector coefficient vector
    f -- function to be interpolated
    refnums -- int number of refinement steps

    Return nothing
    """
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

## @section py_combigrid_grid_converter_2 Regular sparse grids to
## regular combination technique and back: In this example we define a
## regular sparse grid function with a piecewise $d$-linear hat basis
## that interpolates the normal parabola. Then, we transform the
## sparse grid to the corresponding anisotropic grids in the
## combination technique and vice versa. We evaluate the resulting
## surrogates at n randomly chosen samples on the unit hypercube and
## make sure that all of the surrogates we obtain by conversion are
## equal.
def regularGridToRegularGrid(numDims,
                             level,
                             f,
                             n=1000,
                             plot=False,
                             verbose=False):
    """
    Converts a regular sparse grid function to a sparse grid in the
    combination technique and back.

    Arguments:
    numDims -- int number of dimensions
    level -- level of the sparse grid
    f -- function to be interpolated
    n -- int number of random samples on which we evaluate the different sparse grid
         functions to validate the grid conversion
    plot -- bool whether the sparse grid functions are plotted or not (just for numDims=1)
    verbose -- bool verbosity
    """
    ## We generate a iid of uniform samples, which we are going to use to
    ## validate the grid conversion
    x = np.random.rand(n, numDims)
    parameters = DataMatrix(x)

    ## We create a regular sparse grid as usual...
    grid = Grid.createLinearGrid(numDims)
    grid.getGenerator().regular(level)
    alpha = interpolate(grid, f)

    ## We apply now both methods of the grid conversion.
    treeStorage_all = convertHierarchicalSparseGridToCombigrid(grid.getStorage(),
                                                               GridConversionTypes_ALLSUBSPACES)
    treeStorage_full = convertHierarchicalSparseGridToCombigrid(grid.getStorage(),
                                                                GridConversionTypes_COMPLETESUBSPACES)

    ## Note, that we do the conversion just based on the grid points. With
    ## this approach you can easily use different basis function types on
    ## the same grid. We initialize the CombigridOperation on a grid that
    ## spans the same function space as the original hierarchical sparse
    ## grid: hat basis on an equidistant grids without boundary points.
    func = multiFunc(f)
    opt_full = CombigridMultiOperation.createExpUniformLinearInterpolation(numDims, func)
    opt_all = CombigridMultiOperation.createExpUniformLinearInterpolation(numDims, func)

    ## The CombigridOperation expects the points at which you want to
    ## evaluate the interpolant as DataMatrix with the shape (numDims
    ## x numSamples). We, therefore, need to transpose the samples and
    ## initialize the multi operation with them. To set the level
    ## structure we initialize the level manager of the operation with
    ## the storage we have obtained after the conversion.
    parameters.transpose()
    opt_full.setParameters(parameters)
    opt_full.getLevelManager().addLevelsFromStructure(treeStorage_full)
    opt_all.setParameters(parameters)
    opt_all.getLevelManager().addLevelsFromStructure(treeStorage_all)
    parameters.transpose()

    ## If you want you can examine the levels of the combination
    ## technique...
    if verbose:
        print "-" * 80
        print "just full levels:"
        print opt_full.getLevelManager().getSerializedLevelStructure()
        print "-" * 80
        print "all levels:"
        print opt_all.getLevelManager().getSerializedLevelStructure()
        print "-" * 80

    ## We start to transform the grids from the combination technique
    ## back to their hierarchical formulation. We, again, create a
    ## grid with a piecewise $d$-linear basis and initialize the grid
    ## points in its storage by the ones available in the levels of
    ## the combination technique. We do it first for the combination
    ## grids that just contain just those levels where the original
    ## sparse grid had complete subpsaces...
    grid_full = Grid.createLinearGrid(numDims)
    treeStorage_full = opt_full.getLevelManager().getLevelStructure()
    convertCombigridToHierarchicalSparseGrid(treeStorage_full, grid_full.getStorage())

    ## ... and do the same for the version where we considered all
    ## subspaces where at least one grid point was located.
    grid_all = Grid.createLinearGrid(numDims)
    treeStorage_all = opt_all.getLevelManager().getLevelStructure()
    convertCombigridToHierarchicalSparseGrid(treeStorage_all, grid_all.getStorage())

    alpha_full = interpolate(grid_full, f)
    alpha_all = interpolate(grid_all, f)

    ## evaluate the surrogate functions
    y_sg_regular = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid, parameters).eval(alpha, y_sg_regular)

    y_sg_all = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid_all, parameters).eval(alpha_all, y_sg_all)

    y_sg_full = DataVector(parameters.getNrows())
    createOperationMultipleEval(grid_full, parameters).eval(alpha_full, y_sg_full)

    y_ct_all = opt_all.getResult()
    y_ct_full = opt_full.getResult()

    ## tests
    y_sg_regular = y_sg_regular.array().flatten()
    y_ct_all = y_ct_all.array().flatten()
    y_ct_full = y_ct_full.array().flatten()
    y_sg_all = y_sg_all.array().flatten()
    y_sg_full = y_sg_full.array().flatten()

    if plot and numDims == 1:
        x = x.flatten()
        ixs = np.argsort(x)
        plt.figure()
        plt.plot(x[ixs], y_sg_regular[ixs], label="sg regular")
        plt.plot(x[ixs], y_sg_all[ixs], label="sg all")
        plt.plot(x[ixs], y_ct_full[ixs], label="ct full")
        plt.plot(x[ixs], y_ct_all[ixs], label="ct all")
        plt.legend()
        plt.show()

    ## all the functions should be equivalent
    assert np.sum((y_ct_full - y_ct_all) ** 2) < 1e-14
    assert np.sum((y_ct_full - y_sg_regular) ** 2) < 1e-14
    assert np.sum((y_sg_regular - y_sg_all) ** 2) < 1e-14
    assert np.sum((y_sg_regular - y_sg_full) ** 2) < 1e-14
    assert grid_full.getSize() == grid.getSize()
    assert grid_all.getSize() == grid.getSize()


def adaptiveGridToRegularGrid(numDims, level, refnums, f, n=1000, plot=False, verbose=False):
    ## We generate a iid of uniform samples, which we are going to use to
    ## validate the grid conversion
    x = np.random.rand(n, numDims)
    parameters = DataMatrix(x)

    ## We create a regular sparse grid as usual...
    grid = Grid.createLinearGrid(numDims)
    grid.getGenerator().regular(level)
    alpha = interpolate(grid, f)

    grid_adaptive = grid.clone()
    alpha_adaptive = DataVector(alpha)
    refineGrid(grid_adaptive, alpha_adaptive, f, refnums)

    ## We apply now both methods of the grid conversion
    treeStorage_all = convertHierarchicalSparseGridToCombigrid(grid_adaptive.getStorage(),
                                                               GridConversionTypes_ALLSUBSPACES)
    treeStorage_full = convertHierarchicalSparseGridToCombigrid(grid_adaptive.getStorage(),
                                                                GridConversionTypes_COMPLETESUBSPACES)

    ## Note, that we do the conversion just based on the grid points. With
    ## this approach you can easily use different basis function types on
    ## the same grid. We initialize the CombigridOperation on a grid that
    ## spans the same function space as the original hierarchical sparse
    ## grid: hat basis on an equidistant grids without boundary points.
    func = multiFunc(f)
    opt_full = CombigridMultiOperation.createExpUniformLinearInterpolation(numDims, func)
    opt_all = CombigridMultiOperation.createExpUniformLinearInterpolation(numDims, func)

    ## The CombigridOperation expects the points at which you want to
    ## evaluate the interpolant as DataMatrix with the shape (numDims x
    ## numSamples). We initialize the levels for both operations.

    parameters.transpose()
    opt_full.setParameters(parameters)
    opt_full.getLevelManager().addLevelsFromStructure(treeStorage_full)
    opt_all.setParameters(parameters)
    opt_all.getLevelManager().addLevelsFromStructure(treeStorage_all)
    parameters.transpose()

    if verbose:
        print "-" * 80
        print "just full levels:"
        print opt_full.getLevelManager().getSerializedLevelStructure()
        print "-" * 80
        print "all levels:"
        print opt_all.getLevelManager().getSerializedLevelStructure()
        print "-" * 80

    grid_full = Grid.createLinearGrid(numDims)
    treeStorage_full = opt_full.getLevelManager().getLevelStructure()
    convertCombigridToHierarchicalSparseGrid(treeStorage_full, grid_full.getStorage())

    grid_all = Grid.createLinearGrid(numDims)
    treeStorage_all = opt_all.getLevelManager().getLevelStructure()
    convertCombigridToHierarchicalSparseGrid(treeStorage_all, grid_all.getStorage())

    alpha_full = interpolate(grid_full, f)
    alpha_all = interpolate(grid_all, f)

    ## evaluate the surrogate functions
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

    ## tests
    y_sg_regular = y_sg_regular.array().flatten()
    y_sg_adaptive = y_sg_adaptive.array().flatten()
    y_ct_all = y_ct_all.array().flatten()
    y_ct_full = y_ct_full.array().flatten()
    y_sg_all = y_sg_all.array().flatten()
    y_sg_full = y_sg_full.array().flatten()

    if plot and numDims == 1:
        x = x.flatten()
        ixs = np.argsort(x)
        plt.figure()
        plt.plot(x[ixs], y_sg_regular[ixs], label="sg regular")
        plt.plot(x[ixs], y_sg_adaptive[ixs], label="sg adaptive")
        plt.plot(x[ixs], y_ct_full[ixs], label="ct full")
        plt.plot(x[ixs], y_ct_all[ixs], label="ct all")
        plt.plot(x[ixs], y_sg_full[ixs], label="sg full")
        plt.plot(x[ixs], y_sg_all[ixs], label="sg all")
        plt.legend()
        plt.show()

    ## asserts
    assert grid_adaptive.getSize() > grid.getSize()
    assert grid_full.getSize() <= grid_adaptive.getSize()
    assert grid_all.getSize() >= grid.getSize()

    if grid_full.getSize() < grid_all.getSize():
        assert np.sum((y_ct_full - y_ct_all) ** 2) > 1e-14
        assert np.sum((y_sg_regular - y_sg_all) ** 2) > 1e-14

    if grid_full.getSize() == grid.getSize():
        assert np.sum((y_ct_full - y_sg_regular) ** 2) < 1e-14
        assert np.sum((y_sg_full - y_sg_regular) ** 2) < 1e-14


if __name__ == '__main__':
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--numDims', default=4, type=int, help='number of dimensions')
    parser.add_argument('--level', default=4, type=int, help='sparse grid level')
    parser.add_argument('--refnums', default=0, type=int, help='number of refinement steps')
    parser.add_argument('--plot', default=False, action='store_true', help='plot stuff')
    parser.add_argument('--verbose', default=False, action='store_true', help='verbosity')
    args = parser.parse_args()

    if args.refnums == 0:
        regularGridToRegularGrid(args.numDims,
                                 args.level,
                                 f,
                                 plot=args.plot,
                                 verbose=args.verbose)
    else:
        adaptiveGridToRegularGrid(args.numDims,
                                  args.level,
                                  args.refnums,
                                  f,
                                  plot=args.plot,
                                  verbose=args.verbose)
