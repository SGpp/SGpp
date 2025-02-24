#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_combigrid_py Combigrid Example (Python)
##
## In this example, we use the combigrid module to interpolate a test function on a two-dimensional
## regular sparse grid with the combination technique and hierarchical B-splines.
##
## First, we import the required modules.
import pysgpp
import numpy as np
import sys

# skip plotting if Matplotlib cannot be imported (e.g., not installed or no GUI available)
try:
  import matplotlib.pyplot as plt
  from mpl_toolkits.mplot3d import Axes3D
  doPlot = True
except ImportError:
  doPlot = False

## Next, we define a helper function for plotting the resulting functions.
def plotFunction(opEval, surpluses, X):
  if not doPlot: return

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
  ax = fig.add_subplot(projection="3d")
  ax.plot_surface(XX0, XX1, YY)
  ax.plot(X[:,0], X[:,1], "k.", zs=f(X[:,0], X[:,1]), ms=10)

## We define some parameters such as dimensionality and level of the regular sparse grid.
# dimensionality
dim = 2
# regular level
n = 4
# B-spline degree
p = 3
# whether there are points on the boundary
hasBoundary = True
# test function
f = lambda XX0, XX1: np.sin(7*XX0-3)*np.cos(5*XX1-5)

# disable log output
pysgpp.Printer.getInstance().setVerbosity(-1)

## The basis functions are defined via an sgpp::combigrid::HeterogeneousBasis object. In contrast
## to sgpp::base::Basis, this allows for different types of basis functions for the different
## dimensions. However, for this example, we do not need this flexibility, so we use the same
## basis function types for both dimensions.
basis1d = pysgpp.SBsplineBase(p)
basis = pysgpp.HeterogeneousBasis(dim, basis1d)

## An sgpp::combigrid::CombinationGrid is a collection of full grids (nodal subspaces) together
## with scalar-valued coefficients. Here, we construct an sgpp::combigrid::CombinationGrid object
## for a regular sparse grid via the combination technique.
combiGrid = pysgpp.CombinationGrid.fromRegularSparse(dim, n, basis, hasBoundary)

## We obtain the grid points of the regular sparse grid by combining the grid points of all
## full grids that are contained in the combination grid.
gridStorage = pysgpp.HashGridStorage(dim)
combiGrid.combinePoints(gridStorage)

# convert grid storage to array with coordinates of grid points
X = np.array([[gridStorage.getPoint(k).getStandardCoordinate(d) for d in range(dim)]
              for k in range(gridStorage.getSize())])

# evaluate test function at grid points
fX = pysgpp.DataVector(f(X[:,0], X[:,1]))

## We now want to perform an operation on each full grid. For this, we distribute the values of
## the combined grid (sparse grid) to the full grids. The result is a \c std::vector of
## sgpp::base::DataVector; each \c DataVector contains the values at all grid points for one
## specific full grid.
values = pysgpp.DataVectorVector()
combiGrid.distributeValuesToFullGrids(gridStorage, fX, values)

## The operation we want to perform on each full grid is hierarchization. Since the grids are
## full grids, we can use the unidirectional principle for this, which performs 1D hierarchization
## on each pole (one-dimensional sub-grid), iterating over all dimensions.
# copy the values (surpluses will be modified in-place)
surpluses = pysgpp.DataVectorVector(values)

# create pole operation
opPole = pysgpp.OperationPoleVector()
pysgpp.OperationPoleHierarchisationGeneral.fromHeterogenerousBasis(basis, opPole)

# create operation for unidirectional principle and hierarchize in-place
opHier = pysgpp.OperationUPCombinationGrid(combiGrid, opPole)
opHier.apply(surpluses)

## The resulting surpluses are also a \c std::vector of sgpp::base::DataVector, separated by full
## grids. We could combine the full grid surpluses via the combination formula to the sparse grid
## surpluses via \c combineSparseGridValues. However, the operation
## sgpp::combigrid::OperationEvalCombinationGrid does this automatically.
##
## We evaluate the combined function (combination of all full grid interpolants) at some arbitrary
## point, print the value, and plot the function.
# test point at which to evaluate
x = [0.12, 0.34]
xDv = pysgpp.DataVector(x)
print("Value of test function at {}: {:.6g}".format(np.array(x), f(*x)))
# create operation for evaluating and evaluate
opEval = pysgpp.OperationEvalCombinationGrid(combiGrid)
y = opEval.eval(surpluses, xDv)
print("Value of combined sparse grid interpolant at {}: {:.6g}".format(np.array(x), y))
# plot
plotFunction(opEval, surpluses, X)

## Finally, we do the same for one full grid of the combination grid: We evaluate and plot the
## corresponding interpolant. We extract the surpluses from the already calculated \c vector of
## \c DataVector. Alternatively, we could also apply sgpp::combigrid::OperationUPFullGrid with
## opPole to obtain the surpluses for this single full grid.
# select the second full grid of the combination grid (arbitrary choice)
fullGridIndex = 1
fullGrid = combiGrid.getFullGrids()[fullGridIndex]
l = fullGrid.getLevel()
print("Level of selected full grid with index {}: {}".format(fullGridIndex, np.array(l)))

# create operation for evaluating and evaluate
opEval = pysgpp.OperationEvalFullGrid(fullGrid)
y = opEval.eval(surpluses[fullGridIndex], xDv)
print("Value of full grid interpolant at {}: {:.6g}".format(np.array(x), y))

# compute grid points of full grid
X = pysgpp.DataMatrix(0, 0)
pysgpp.IndexVectorRange.getPoints(fullGrid, X)

# convert resulting sgpp::base::DataMatrix to NumPy array
X = np.array([[X.get(k, j) for j in range(X.getNcols())] for k in range(X.getNrows())])

# plot
plotFunction(opEval, surpluses[fullGridIndex], X)

if doPlot: plt.show()
else: print("Skipping plots due to failed import of Matplotlib.")

## The example program outputs the following results:
## \verbinclude combigrid.output.txt
##
## We see that the value of the combined sparse grid interpolant at the evaluation point is
## closer to the actual value of the test function than the value of the chosen full grid
## interpolant, which corresponds to the full grid of level \f$(3, 1)\f$.
