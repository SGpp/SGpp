# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import binom

from pysgpp import Grid, DataVector, SurplusRefinementFunctor, createOperationHierarchisation, RegularizationType_Laplace
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d, plotSG2d, \
    plotGrid2d
from pysgpp.extensions.datadriven.uq.dists import MultivariateNormal, SGDEdist, Beta, J
from pysgpp.extensions.datadriven.uq.operations import (hierarchize,
                                                        evalSGFunctionMulti)
from pysgpp.extensions.datadriven.uq.operations.forcePositivity import OperationMakePositiveFast, InterpolateFunction, ScaledMinOfParents
from pysgpp.extensions.datadriven.uq.quadrature import doQuadrature
from pysgpp.extensions.datadriven.uq.dists.Lognormal import Lognormal
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.findIntersectionsSubspaceBased import IntersectionSubspaceCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.findIntersections import IntersectionCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.fullGridSearch import FullGridCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.searchNextLevel import SearchLevelWiseForCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.localFullGridSearch import LocalFullGridCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.localHierarchicalIntersectionSearch import LocalHierarchicalIntersectionCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.estimateDensity import EstimateDensityAlgorithm
from pysgpp import createOperationMakePositive, \
    MakePositiveCandidateSearchAlgorithm_Intersections, \
    MakePositiveInterpolationAlgorithm_SetToZero, \
    MakePositiveInterpolationAlgorithm_InterpolateExp, \
    MakePositiveInterpolationAlgorithm_InterpolateBoundaries1d, SizeList, \
    HashGridPoint, MakePositiveCandidateSearchAlgorithm_HybridFullIntersections, \
    GridType_Linear, GridType_Poly, RegularGridConfiguration
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import checkPositivity
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotSG3d, plotDensity3d
from pysgpp.pysgpp_swig import GridType_PolyBoundary, \
    MakePositiveCandidateSearchAlgorithm_FullGrid, GridType_LinearClenshawCurtis, \
    MakePositiveCandidateSearchAlgorithm_IntersectionsJoin
from pysgpp.extensions.datadriven.uq.dists.Normal import Normal


# parameters
gridConfig = RegularGridConfiguration()
# gridConfig.type_ = GridType_Linear
gridConfig.type_ = GridType_Linear
gridConfig.boundaryLevel_ = 0
# issues with: d = 5, l = 3, ref = 4
numDims = 6
level = 4
refnums = 0
consistentGrid = False
# candidateSearchAlgorithm = MakePositiveCandidateSearchAlgorithm_FullGrid
# candidateSearchAlgorithm = MakePositiveCandidateSearchAlgorithm_Intersections
candidateSearchAlgorithm = MakePositiveCandidateSearchAlgorithm_IntersectionsJoin
# candidateSearchAlgorithm = MakePositiveCandidateSearchAlgorithm_HybridFullIntersections
# interpolationAlgorithm = MakePositiveInterpolationAlgorithm_InterpolateBoundaries1d
interpolationAlgorithm = MakePositiveInterpolationAlgorithm_SetToZero
plot = True
verbose = True
code = "c++"

gridConfig.dim_ = numDims
gridConfig.level_ = level
gridConfig.maxDegree_ = level + 1

mu = np.ones(numDims) * 0.5
cov = np.diag(np.ones(numDims) * 0.1 / 10.)

dist = J([Normal(0.5, 1./ 16., 0, 1)] * numDims)
# dist = MultivariateNormal(mu, cov, 0, 1)  # problems in 3d/l2
# dist = J([Beta(5, 4, 0, 1)] * numDims)  # problems in 5d/l3
# dist = J([Lognormal(0.2, 0.7, 0, 1)] * numDims)  # problems in 5d/l3

trainSamples = dist.rvs(1000)
testSamples = dist.rvs(1000)

# plot analytic density
if numDims == 2 and plot:
    fig = plt.figure()
    plotDensity2d(dist)
    plt.title("analytic, kldivergence = %g" % dist.klDivergence(dist, testSamples))
    fig.show()

    fig, ax = plotDensity3d(dist)
    ax.set_xlabel(r"$x_1$")
    ax.set_ylabel(r"$x_2$")
    ax.set_zlabel(r"$f(x_1, x_2)$")
    fig.show()

# get a sparse grid approximation
grid = Grid.createGrid(gridConfig)
gridGen = grid.getGenerator()
gridGen.regular(level)
gs = grid.getStorage()

# now refine adaptively 5 times
p = DataVector(gs.getDimension())
alpha = DataVector(gs.getSize())

# set function values in alpha
for i in range(gs.getSize()):
    gs.getPoint(i).getStandardCoordinates(p)
    alpha[i] = dist.pdf(p.array())

# hierarchize
createOperationHierarchisation(grid).doHierarchisation(alpha)

for refnum in range(refnums):
    # refine a single grid point each time
    gridGen.refine(SurplusRefinementFunctor(alpha, 1))
    print("refinement step {}, new grid size: {}".format(refnum + 1, gs.getSize()))

    # extend alpha vector (new entries uninitialized)
    alpha.resize(gs.getSize())

    # set function values in alpha
    for i in range(gs.getSize()):
        gs.getPoint(i).getStandardCoordinates(p)
        alpha[i] = dist.pdf(p.array())

    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)


alpha = alpha.array()
sgdeDist = SGDEdist(grid, alpha,
                    trainData=trainSamples, bounds=dist.getBounds())
print("l=%i: (gs=%i) -> %g (%g, %g)," % (level,
                                         sgdeDist.grid.getSize(),
                                         dist.klDivergence(sgdeDist, testSamples),
                                         sgdeDist.crossEntropy(testSamples),
                                         sgdeDist.vol))
print("-" * 80)

if numDims == 2 and plot:
    # plot the result
    fig = plt.figure()
    plotGrid2d(grid, alpha, show_numbers=False)
#     plt.title("neg: #gp = %i, kldivergence = %g, log = %g" % (grid.getStorage().getSize(),
#                                                               dist.klDivergence(sgdeDist, testSamples),
#                                                               dist.crossEntropy(testSamples)))
    fig.show()

    fig, ax, _ = plotSG3d(grid, sgdeDist.alpha)
    ax.set_title("negative")
    fig.show()

C = 0
M = np.sum([1 for i in range(len(alpha)) if alpha[i] < 0])
for d in range(2, numDims + 1):
    C += binom(M, d)
print("predicted comparison costs = %i" % C)
print("full grid                  = %i" % ((2 ** level - 1) ** numDims,))

if code == "c++":
    alpha_vec = DataVector(alpha)
    opMakePositive = createOperationMakePositive(candidateSearchAlgorithm,
                                                 interpolationAlgorithm,
                                                 consistentGrid, verbose)
    opMakePositive.makePositive(grid, alpha_vec)
    alpha = alpha_vec.array()
elif code == "python":
    cand = None
    cand = IntersectionCandidates()
#     cand = IntersectionSubspaceCandidates()
    # cand = FullGridCandidates(grid)
    # cand = LocalFullGridCandidates(grid)
    # cand = SearchLevelWiseForCandidates(grid, alpha)
    # cand = LocalHierarchicalIntersectionCandidates(grid)

    alg = None
    # alg = ScaledMinOfParents()
    # alg = EstimateDensityAlgorithm(trainSamples)

    # now make the function positive
    opPositive = OperationMakePositiveFast(grid,
                                           candidateSearchAlgorithm=cand,
                                           interpolationAlgorithm=alg)

    grid, alpha = opPositive.makePositive(alpha)
else:
    raise AttributeError("code '%s' not known" % code)

# # check for minimality
# addedGridPoints = opMakePositive.getAddedGridPointsForPositivity()
# x = DataVector(numDims)
# for ix in addedGridPoints:
#     coarsedGrid = grid.clone()
#     coarsedGridStorage = coarsedGrid.getStorage()
#     gp = HashGridPoint(coarsedGridStorage.getPoint(ix))
#     toBeRemoved = ListIndex()
#     toBeRemoved.push_back(ix)
#     remainingIndex = coarsedGridStorage.deletePoints(toBeRemoved)
#     coarsedAlpha = DataVector(alpha)
#     coarsedAlpha.restructure(remainingIndex)
#
#     gp.getStandardCoordinates(x);
#     opEvalNaive = createOperationEvalNaive(*coarsedGrid)
#     value = opEvalNaive.eval(coarsedAlpha, x)
#
#     print value

# security check for positiveness
neg = checkPositivity(grid, alpha)

if len(neg) > 0:
    print("warning: the sparse grid function is not positive")
#             raise AttributeError("the sparse grid function is not positive")
    # check at which grid points the function is negative
#             for i, (yi, gp) in neg.items():
#                     print "|%s|_1 = %i, %s -> %g" % ([gp.getLevel(d) for d in xrange(numDims)],
#                                                      np.sum([gp.getLevel(d) for d in xrange(numDims)]),
#                                                      [gp.getIndex(d) for d in xrange(numDims)],
#                                                      yi)

sgdeDist = SGDEdist(grid, alpha,
                    trainData=trainSamples, bounds=dist.getBounds())
print("-" * 80)
print("(gs=%i) -> %g (%g, %g)" % (sgdeDist.grid.getSize(),
                                  dist.klDivergence(sgdeDist, testSamples),
                                  sgdeDist.crossEntropy(testSamples),
                                  sgdeDist.vol))

if numDims == 2 and plot:
    # plot the result
#     fig = plt.figure()
#     plotGrid2d(grid, alpha, show_numbers=False)
#     plt.title("pos: #gp = %i, kldivergence = %g, log = %g" % (grid.getStorage().getSize(),
#                                                               dist.klDivergence(sgdeDist, testSamples),
#                                                               dist.crossEntropy(testSamples)))
    fig.show()

    fig = plt.figure()
    plotSG2d(grid, sgdeDist.alpha, show_negative=True, show_grid_points=True)
    fig.show()
    fig, ax, _ = plotSG3d(grid, sgdeDist.alpha)
#     ax.set_title("positive")
    fig.show()
        
    plt.show()

