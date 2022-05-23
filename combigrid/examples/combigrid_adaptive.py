#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_combigrid_adaptive_py Combigrid Example Dimensional Adaptivity (Python)
##
## In this example, we use the combigrid module for dimensional adaptivity
##
##
## First, we import the required modules.
import pysgpp
import numpy as np
import sys

# skip plotting if Matplotlib cannot be imported (e.g., not installed or no GUI available)
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from pysgpp.extensions.combigrid.plotDeltas3d import plotDeltas3D
    doPlot = True
except ImportError:
    doPlot = False


# we define a helper function for plotting the resulting functions.
def plotFunction(opEval, surpluses, X):
    if not doPlot:
        return

    # generate a meshgrid for plotting
    xx0 = np.linspace(0, 1, 65)
    xx1 = np.linspace(0, 1, 65)
    XX0, XX1 = np.meshgrid(xx0, xx1)
    XX = pysgpp.DataMatrix(np.column_stack([XX0.flatten(), XX1.flatten()]))

    # evaluate interpolant at meshgrid
    YY = pysgpp.DataVector(0)
    opEval.multiEval(surpluses, XX, YY)

    # convert resulting sgpp::base::DataVector to NumPy array
    YY = np.reshape(np.array([YY[k] for k in range(YY.getSize())]), XX0.shape)

    # actual plotting
    fig = plt.figure(figsize=(6, 6))
    ax = fig.gca(projection="3d")
    ax.plot_surface(XX0, XX1, YY)
    ax.plot(X[:, 0], X[:, 1], "k.", zs=f(X[:, 0], X[:, 1]), ms=10)


# We define parameters and perform hierarchization just as in the combigrid example
dim = 2
n = 4
p = 3
hasBoundary = True
# test function


def f(XX0, XX1): return np.sin(7*XX0-3)*np.cos(5*XX1-5)


# disable log output
pysgpp.Printer.getInstance().setVerbosity(-1)

basis1d = pysgpp.SBsplineBase(p)
basis = pysgpp.HeterogeneousBasis(dim, basis1d)
combiGrid = pysgpp.CombinationGrid.fromRegularSparse(
    dim, n, basis, hasBoundary)

gridStorage = pysgpp.HashGridStorage(dim)
combiGrid.combinePoints(gridStorage)

X = np.array([[gridStorage.getPoint(k).getStandardCoordinate(d) for d in range(dim)]
              for k in range(gridStorage.getSize())])

fX = pysgpp.DataVector(f(X[:, 0], X[:, 1]))

values = pysgpp.DataVectorVector()
combiGrid.distributeValuesToFullGrids(gridStorage, fX, values)

surpluses = pysgpp.DataVectorVector(values)
opPole = pysgpp.OperationPoleVector()
pysgpp.OperationPoleHierarchisationGeneral.fromHeterogenerousBasis(
    basis, opPole)
opHier = pysgpp.OperationUPCombinationGrid(combiGrid, opPole)
opHier.apply(surpluses)


# Let's assume that this time we are even more interested in the interpolated value
# at a given point x, so much so that we would like to extend the combination scheme
# in a way that makes this value more accurate.
# f(x) can be considered the quantity of interest or QoI.

x = [0.12, 0.34]
xDv = pysgpp.DataVector(x)
print("Value of test function at {}: {:.6g}".format(np.array(x), f(*x)))
# create operation for evaluating and evaluate
opEval = pysgpp.OperationEvalCombinationGrid(combiGrid)
y = opEval.eval(surpluses, xDv)
print("Value of combined sparse grid interpolant at {}: {:.6g}".format(np.array(x), y))
# plot
plotFunction(opEval, surpluses, X)


# we generate a relevance calculator (for known values) and a priority estimator (for unknown values)
# and test them
weightedRelevanceCalculator = pysgpp.WeightedRelevanceCalculator()
relevance = weightedRelevanceCalculator.calculate(
    pysgpp.LevelVector([4, 5, 5]), 0.8)


averagingPriorityEstimator = pysgpp.AveragingPriorityEstimator()
priority = averagingPriorityEstimator.estimatePriority(pysgpp.LevelVector([5, 5, 5]),
                                                       pysgpp.map_levelvector_real({
                                                           pysgpp.LevelVector([5, 4, 5]): 0.8,
                                                           pysgpp.LevelVector([4, 5, 5]): 0.8
                                                       }))

# for each grid in our combination scheme, we calculate the value at point x
subgridValuesAtX = pysgpp.DoubleVector(len(combiGrid.getFullGrids()), 0.)
for fullGridIndex in range(len(combiGrid.getFullGrids())):
    fullGrid = combiGrid.getFullGrids()[fullGridIndex]
    l = fullGrid.getLevel()
    print("Level of selected full grid with index {}: {}".format(
        fullGridIndex, np.array(l)))

    # create operation for evaluating and evaluate
    opEval = pysgpp.OperationEvalFullGrid(fullGrid)
    y = opEval.eval(surpluses[fullGridIndex], xDv)
    print("Value of full grid interpolant at {}: {:.6g}".format(np.array(x), y))
    subgridValuesAtX[fullGridIndex] = y

# we start an AdaptiveCombinationGridGenerator from the full grids that are in the combiGrid we already have
adaptiveCombinationGridGenerator = pysgpp.AdaptiveCombinationGridGenerator.fromCombinationGrid(
    combiGrid, subgridValuesAtX)
# same as:
# adaptiveCombinationGridGenerator = pysgpp.AdaptiveCombinationGridGenerator.fromCombinationGrid(combiGrid, subgridValuesAtX,
#                                                                                               weightedRelevanceCalculator, averagingPriorityEstimator)

# all of the full grids so far should already be in the generator's 'old set', so this call here does not change anything:
adaptiveCombinationGridGenerator.adaptAllKnown()
oldSet = adaptiveCombinationGridGenerator.getOldSet()

# look at the current result of the combination -- same as the combigrid evaluation above
adaptiveCombinationGridGenerator.getCurrentResult()

levels = pysgpp.LevelVectorVector()
for fullGridIndex in range(len(combiGrid.getFullGrids())):
    fullGrid = combiGrid.getFullGrids()[fullGridIndex]
    l = fullGrid.getLevel()
    levels.push_back(l)
adaptiveCombinationGridGenerator.getDeltas(levels)


# looking at the 'active set', we see which full grid spaces could potentially be interesting next,
# they are the upper neighbors of the old set
activeSet = adaptiveCombinationGridGenerator.getActiveSet()

# some of them are more interesting than others
priorities = adaptiveCombinationGridGenerator.getPriorities()

# we ask what the next most sensible subspaces to evaluate should be
queue = adaptiveCombinationGridGenerator.getPriorityQueue()

# we calculate the value at x in the corresponding full grids
newValue = {}
for qu in queue:
    q = pysgpp.LevelVector(qu)

    # create full grid, interpolate and hierarchize
    fullGrid = pysgpp.FullGrid(q, basis)
    gridRange = pysgpp.IndexVectorRange(fullGrid)
    value = pysgpp.DataVector(fullGrid.getNumberOfIndexVectors())

    rangeIterator = gridRange.begin()
    while rangeIterator != gridRange.end():
        #         print(*rangeIterator.__deref__())
        z = pysgpp.DataVector(dim)
        rangeIterator.getStandardCoordinates(z)
        value[rangeIterator.getSequenceNumber()] = f(z[0], z[1])
        try:
            rangeIterator.__next__()
        except StopIteration:
            break

    # store the interpolated values for later reuse
    newValue[tuple(q)] = value
    surplus = pysgpp.DataVector(value)
    opHier = pysgpp.OperationUPFullGrid(fullGrid, opPole)
    opHier.apply(surplus)

    # create operation for evaluating and evaluate
    opEval = pysgpp.OperationEvalFullGrid(fullGrid)
    y = opEval.eval(surplus, xDv)
    print("Value of full grid interpolant at {}: {:.6g}".format(np.array(x), y))

    # and add the results to the generator
    adaptiveCombinationGridGenerator.setQoIInformation(q, y)

# we can test how well the priority estimator had performed
# by looking at the relevance of the now known results
relevances = adaptiveCombinationGridGenerator.getRelevanceOfActiveSet()
priorities, relevances

# we adapt to all grids with known values
adaptiveCombinationGridGenerator.adaptAllKnown()
currentSet = adaptiveCombinationGridGenerator.getOldSet()

# and get the new total result of the adapted combination
adaptiveCombinationGridGenerator.getCurrentResult()

# now, we can also update our combiGrid by adding the newly calculated subspaces
newCombiGrid = adaptiveCombinationGridGenerator.getCombinationGrid(basis)

numberNewLevels = len(newCombiGrid.getFullGrids())
newLevels = pysgpp.LevelVectorVector()
for fullGridIndex in range(numberNewLevels):
    fullGrid = newCombiGrid.getFullGrids()[fullGridIndex]
    l = fullGrid.getLevel()
    newLevels.push_back(l)

newCoefficients = pysgpp.getStandardCoefficientsFromLevelSet(newLevels)

# and by calculating the new surplusses...
newValues = pysgpp.DataVectorVector()
newValues.reserve(numberNewLevels)
for newLevel in newLevels:
    if newLevel in levels:
        # if fg was already there in the old combiGrid, re-use the values
        #         i = levels.index(fg) # this does not currently work, manual "hack" below
        i = 0
        while levels[i] != newLevel and i < numberNewLevels:
            i = i+1
        newValues.push_back(values[i])
    else:
        # use the newly interpolated values
        newValues.push_back(newValue[newLevel])

newSurpluses = pysgpp.DataVectorVector(newValues)
newOpHier = pysgpp.OperationUPCombinationGrid(newCombiGrid, opPole)
newOpHier.apply(newSurpluses)

# we can see that the newly interpolated sparse grid value is in fact
# closer to the analytical solution
opEval = pysgpp.OperationEvalCombinationGrid(newCombiGrid)
y = opEval.eval(newSurpluses, xDv)
print("Value of combined sparse grid interpolant at {}: {:.6g}".format(np.array(x), y))

# and looking at the delta increments contributed by each of the grids,
# we can see which grid resolutions had the biggest impact:
if doPlot:
    plotDeltas3D(adaptiveCombinationGridGenerator, ['l_0', 'l_1'])
    # when you change the dimension in the example to 3 or more, you can see the 3-d visualization of the deltas:
    # plotDeltas3d.plotDeltas3D(adaptiveCombinationGridGenerator)
