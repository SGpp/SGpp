from builtins import range
import numpy as np
import matplotlib.pyplot as plt

from pysgpp import (Grid, DataVector,
                    createOperationLimitFunctionValueRange,
                    MakePositiveCandidateSearchAlgorithm_Intersections,
                    MakePositiveInterpolationAlgorithm_SetToZero)
from pysgpp.extensions.datadriven.uq.plot import plotDensity2d, plotSG2d, \
    plotFunction2d, plotFunction1d
from pysgpp.extensions.datadriven.uq.operations import (hierarchize,
                                                        evalSGFunctionMulti)
from pysgpp.extensions.datadriven.uq.operations.forcePositivity import OperationMakePositiveFast, InterpolateFunction
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.findIntersections import IntersectionCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.fullGridSearch import FullGridCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.searchNextLevel import SearchLevelWiseForCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.localFullGridSearch import LocalFullGridCandidates
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSG1d
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import addConst
from pysgpp.pysgpp_swig import createOperationMakePositive

numDims = 4
level = 4
code = "python"


def hierarchizeFun(fun):
    nodalValues = np.ndarray(grid.getSize())
    p = DataVector(gs.getDimension())
    for i in range(gs.getSize()):
        gs.getPoint(i).getStandardCoordinates(p)
        nodalValues[i] = fun(p.array())

    return hierarchize(grid, nodalValues)

def sin(x):
    return np.sum(np.sin(np.pi * (2 * x - 1.0))) + 0.8 * len(x)

fun = sin

# plot analytic function
if numDims == 1:
    fig = plt.figure()
    plotFunction1d(fun)
    plt.title("analytic")
    fig.show()
elif numDims == 2:
    fig = plt.figure()
    plotFunction2d(fun)
    plt.title("analytic")
    fig.show()

# get a sparse grid approximation
grid = Grid.createLinearGrid(numDims)
grid.getGenerator().regular(level)
gs = grid.getStorage()

alpha = hierarchizeFun(fun)
print("l=%i: (gs=%i)" % (level, grid.getSize()))
print("-" * 80)

# plot the result
if numDims < 3:
    fig = plt.figure()
    if numDims == 1:
        plotSG1d(grid, alpha)
    elif numDims == 2:
        plotSG2d(grid, alpha, show_negative=False, show_grid_points=True)

    plt.title("neg: #gp = %i" % grid.getStorage().getSize())
    fig.show()

if code == "c++":
    newGrid = Grid.createLinearGrid(numDims)
    alpha_vec = DataVector(alpha)
    opLimit = createOperationMakePositive(grid, MakePositiveCandidateSearchAlgorithm_Intersections, MakePositiveInterpolationAlgorithm_SetToZero, True)
    opLimit.makePositive(newGrid, alpha_vec)

    grid = newGrid
    alpha = alpha_vec.array()
else:
    grid, alpha = OperationMakePositiveFast(grid).makePositive(alpha)

    print("-" * 80)
    print("l=%i: (gs=%i)" % (level, grid.getSize()))

# plot the final result
if numDims < 3:
    fig = plt.figure()
    if numDims == 1:
        plotSG1d(grid, alpha)
    elif numDims == 2:
        plotSG2d(grid, alpha, show_negative=False, show_grid_points=True)

    plt.title("pos: #gp = %i" % grid.getStorage().getSize())
    fig.show()

    plt.show()
