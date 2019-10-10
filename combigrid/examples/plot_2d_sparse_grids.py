#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_plot_2d_sparse_grids_py plot_2d_sparse_grids.py
## plots anisotropic full grids that form part of the combination technique

try:
    from argparse import ArgumentParser
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    import pysgpp
    from pysgpp.pysgpp_swig import DataVector, CombigridOperation
    from pysgpp.extensions.datadriven.uq.dists import J, Beta, Uniform, Normal
    from pysgpp.extensions.datadriven.uq.plot.colors import load_color, savefig,\
        load_font_properties
    from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d

except ImportError as e:
    print(e.__class__.__name__ + ": " + e.msg)
    print("Skipping example...")
    exit(0)


def g(x):
    return np.prod([4 * xi * (1 - xi) for xi in x.array()])


# We have to wrap f in a pysgpp.MultiFunction object.
func = pysgpp.multiFunc(g)
numDims = 2

if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('--level', default=2, type=int, help="minimum level of regular grids")
    parser.add_argument('--marginalType', default="beta", type=str, help="marginals")
    args = parser.parse_args()

    if args.marginalType == "uniform":
        marginal = Uniform(0, 1)
    elif args.marginalType == "beta":
        marginal = Beta(5, 10)
    else:
        marginal = Normal(0.5, 0.1, 0, 1)

    # plot pdf
    dist = J([marginal] * numDims)
    fig = plt.figure()
    plotDensity2d(dist)
    savefig(fig, "/tmp/%s" % (args.marginalType,))
    plt.close(fig)

    w = pysgpp.singleFunc(marginal.pdf)

    grids = pysgpp.AbstractPointHierarchyVector()
    grids.push_back(pysgpp.CombiHierarchies.linearLeja(w))
    grids.push_back(pysgpp.CombiHierarchies.linearLeja(w))

    evaluators = pysgpp.FloatScalarAbstractLinearEvaluatorVector()
    evaluators.push_back(pysgpp.CombiEvaluators.polynomialInterpolation())
    evaluators.push_back(pysgpp.CombiEvaluators.polynomialInterpolation())

    # To create a CombigridOperation object with our own configuration, we have to provide a
    # LevelManager as well:
    levelManager = pysgpp.WeightedRatioLevelManager()
    operation = pysgpp.CombigridOperation(grids, evaluators, levelManager, func)

    # We can add regular levels like before:
    levelManager.addRegularLevels(args.level)

    # We can also fetch the used grid points and plot the grid:
    grid = levelManager.getGridPointMatrix()
    gridList = [[grid.get(r, c) for c in range(grid.getNcols())] for r in range(grid.getNrows())]

    fig = plt.figure()
    plt.plot(gridList[0], gridList[1], " ",
             color=load_color(0),
             marker='o', markersize=10)
    plt.axis('off')
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0, 0), 1, 1, fill=None, alpha=1, linewidth=2))
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.title(r"Sparse Grid $\ell=%i$" % args.level,
              fontproperties=load_font_properties())
    savefig(fig, "/tmp/sparse_grid_l%i_%s" % (args.level, args.marginalType))
