import numpy as np
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d, plotSG2d
from pysgpp.extensions.datadriven.uq.dists import MultivariateNormal, SGDEdist, Beta, J
from pysgpp import Grid, DataVector
from pysgpp.extensions.datadriven.uq.operations import (hierarchize,
                                                        evalSGFunctionMulti)
from pysgpp.extensions.datadriven.uq.operations.forcePositivity import OperationMakePositiveFast, InterpolateFunction, InterpolateParents
from pysgpp.extensions.datadriven.uq.quadrature import doQuadrature
from pysgpp.extensions.datadriven.uq.dists.Lognormal import Lognormal
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.findIntersections import IntersectionCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.fullGridSearch import FullGridCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.searchNextLevel import SearchLevelWiseForCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.localFullGridSearch import LocalFullGridCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.localHierarchicalIntersectionSearch import LocalHierarchicalIntersectionCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.estimateDensity import EstimateDensityAlgorithm
from pysgpp import RegularizationType_Laplace, SurplusRefinementFunctor, createOperationHierarchisation

# parameters
numDims = 2
level = 2
refnums = 0
plot = True

mu = np.ones(numDims) * 0.5
cov = np.diag(np.ones(numDims) * 0.1 / 5.)

dist = MultivariateNormal(mu, cov, 0, 1)  # problems in 3d/l2
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

# get a sparse grid approximation
sgdeDist = SGDEdist.byLearnerSGDEConfig(trainSamples,
                                        config={"grid_level": level,
                                                "grid_type": "Linear",
                                                "refinement_numSteps": refnums,
                                                "refinement_numPoints": 1,
                                                "regularization_type": "Laplace",
                                                "crossValidation_lambda": 0.000562341,
                                                "crossValidation_enable": True,
                                                "crossValidation_kfold": 5,
                                                "crossValidation_silent": True},
                                        bounds=dist.getBounds())

grid, alpha = sgdeDist .grid, sgdeDist.alpha.array()

print "l=%i: (gs=%i) -> %g (%g, %g)," % (level,
                                         sgdeDist.grid.getSize(),
                                         dist.klDivergence(sgdeDist, testSamples),
                                         sgdeDist.crossEntropy(testSamples),
                                         sgdeDist.vol)
print "-" * 80

if numDims == 2 and plot:
    # plot the result
    fig = plt.figure()
    plotSG2d(grid, sgdeDist.alpha, show_negative=True, show_grid_points=True)
    plt.title("neg: #gp = %i, kldivergence = %g, log = %g" % (grid.getStorage().getSize(),
                                                              dist.klDivergence(sgdeDist, testSamples),
                                                              dist.crossEntropy(testSamples)))
    fig.show()

cand = None
# cand = IntersectionCandidates(grid, alpha)
# cand = FullGridCandidates(grid)
cand = LocalFullGridCandidates(grid)
# cand = SearchLevelWiseForCandidates(grid, alpha)
# cand = LocalHierarchicalIntersectionCandidates(grid)

alg = None
alg = EstimateDensityAlgorithm(trainSamples,
                               lmbd=0.0,  # 0.000562341,
                               regularizationType=RegularizationType_Laplace)
# now make the function positive
opPositive = OperationMakePositiveFast(grid,
                                       candidateSearchAlgorithm=cand,
                                       interpolationAlgorithm=alg)

grid, alpha = opPositive.makePositive(alpha)
# # grid, alpha = makePositiveSkript(grid, alpha, dist.pdf)

sgdeDist = SGDEdist(grid, alpha)
print "-" * 80
print "(gs=%i) -> %g (%g, %g)" % (sgdeDist.grid.getSize(),
                                  dist.klDivergence(sgdeDist, testSamples),
                                  sgdeDist.crossEntropy(testSamples),
                                  sgdeDist.vol)

if numDims == 2 and plot:
    # plot the result
    fig = plt.figure()
    plotSG2d(grid, sgdeDist.alpha, show_negative=True, show_grid_points=True)
    plt.title("pos: #gp = %i, kldivergence = %g, log = %g" % (grid.getStorage().getSize(),
                                                              dist.klDivergence(sgdeDist, testSamples),
                                                              dist.crossEntropy(testSamples)))
    fig.show()
    plt.show()

