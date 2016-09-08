import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib2tikz import save as tikz_save

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
from pysgpp.pysgpp_swig import createOperationMakePositive, \
    MakePositiveCandidateSearchAlgorithm_Intersections, \
    MakePositiveInterpolationAlgorithm_SetToZero, \
    MakePositiveInterpolationAlgorithm_InterpolateExp, \
    MakePositiveInterpolationAlgorithm_InterpolateBoundaries1d, IndexList, \
    HashGridPoint
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import checkPositivity, \
    evalSGFunction, getBoundsOfSupport
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotSG3d, plotDensity3d


def plotIntersections(grid, alpha, (i, gpi), (j, gpj), gpintersection):
    # -----------------------------------------------------------------
    # plot result
    gs = grid.getStorage()

    fig = plt.figure()
    for k in xrange(gs.getSize()):
        gp = gs.getPoint(k)
        x, y = gp.getStandardCoordinate(0), gp.getStandardCoordinate(1)
        if alpha[k] < 0.0:
            plt.plot(x, y, "v ", color="red")
        else:
            plt.plot(x, y, "^ ", color="white")

    # annotate the
    plt.annotate(str(i), (gpi.getStandardCoordinate(0), gpi.getStandardCoordinate(1)))
    plt.annotate(str(j), (gpj.getStandardCoordinate(0), gpj.getStandardCoordinate(1)))

    # draw support
    # get level index
    for gp, col in ((gpi, "b"), (gpj, "r")):
        l0, i0, l1, i1 = gp.getLevel(0), gp.getIndex(0), gp.getLevel(1), gp.getIndex(1)
        xlim0, xlim1 = getBoundsOfSupport(l0, i0), getBoundsOfSupport(l1, i1)
        diff0, diff1 = xlim0[1] - xlim0[0], xlim1[1] - xlim1[0]
        currentAxis = plt.gca()
        currentAxis.add_patch(Rectangle((xlim0[0], xlim1[0]),
                                        diff0, diff1, facecolor=col,
                                        alpha=0.5))

#         xlim = self.computeBoundsOfOverlappingPatch(gpi, gpj)
#         currentAxis = plt.gca()
#         currentAxis.add_patch(Rectangle((xlim[0, 0], xlim[1, 0]),
#                                         np.diff(xlim[0, :2]), np.diff(xlim[1, :2]),
#                                         facecolor="gray", alpha=0.9))

    plt.plot(gpintersection.getStandardCoordinate(0),
             gpintersection.getStandardCoordinate(1),
             "o ", color="black")

    plt.xlim(0, 1)
    plt.ylim(0, 1)

    return fig


# parameters
numDims = 2
level = 2
# interpolationAlgorithm = MakePositiveInterpolationAlgorithm_InterpolateBoundaries1d
interpolationAlgorithm = MakePositiveInterpolationAlgorithm_SetToZero
plot = True
verbose = False

mu = np.ones(numDims) * 0.5
cov = np.diag(np.ones(numDims) * 0.08 / 5.)

dist = MultivariateNormal(mu, cov, 0, 1)  # problems in 3d/l2

trainSamples = dist.rvs(1000)
testSamples = dist.rvs(1000)

# plot analytic density
if numDims == 2 and plot:
    fig, ax = plotDensity3d(dist)
    ax.set_xlabel(r"$x_1$")
    ax.set_ylabel(r"$x_2$")
    ax.set_zlabel(r"$f(x_1, x_2)$")
    fig.show()

# get a sparse grid approximation
grid = Grid.createLinearGrid(numDims)
gridGen = grid.getGenerator()
gridGen.regular(level)
gs = grid.getStorage()

# now refine adaptively 5 times
p = DataVector(gs.getDimension())
alpha = DataVector(gs.getSize())

# set function values in alpha
for i in xrange(gs.getSize()):
    gs.getPoint(i).getStandardCoordinates(p)
    alpha[i] = dist.pdf(p.array())

# hierarchize
createOperationHierarchisation(grid).doHierarchisation(alpha)

alpha = alpha.array()


# create list of missing full grid points
gps = []
gp0 = HashGridPoint(numDims)
gp0.set(0, 2, 1)
gp0.set(1, 2, 1)
gps.append((gp0, 3, gs.getPoint(3), 1, gs.getPoint(1)))

gp1 = HashGridPoint(numDims)
gp1.set(0, 2, 1)
gp1.set(1, 2, 3)
gps.append((gp1, 4, gs.getPoint(4), 1, gs.getPoint(1)))

gp2 = HashGridPoint(numDims)
gp2.set(0, 2, 3)
gp2.set(1, 2, 1)
gps.append((gp2, 3, gs.getPoint(3), 2, gs.getPoint(2)))

gp3 = HashGridPoint(numDims)
gp3.set(0, 2, 3)
gp3.set(1, 2, 3)
gps.append((gp3, 4, gs.getPoint(4), 2, gs.getPoint(2)))

gs = grid.getStorage()
x = DataVector(numDims)

# plot the result
i = 0

# figGrid = plt.figure()
# plotGrid2d(grid, alpha, show_numbers=True)
# plt.title("iteration = %i" % i)
# figGrid.show()

figFun, axFun, _ = plotSG3d(grid, alpha)
axFun.set_title("iteration = %i" % i)
plt.savefig('figures/positivity_normal_i%i.pdf' % i)

for i, (gpintersection, ix, gpi, jx, gpj) in enumerate(gps):
    gs.insert(gpintersection)
    gpintersection.getStandardCoordinates(x)
    alpha = np.append(alpha, -evalSGFunction(grid, alpha, x.array()))

#     # plot the result
#     fig = plotIntersections(grid, alpha, (i, gpi), (j, gpj), gpintersection)
#     fig.show()

    figFun, axFun, _ = plotSG3d(grid, alpha)
    axFun.set_title("iteration = %i" % (i + 1,))
    plt.savefig('figures/positivity_normal_i%i.pdf' % (i + 1,))
    plt.close(fig)

if plot:
    plt.show()
