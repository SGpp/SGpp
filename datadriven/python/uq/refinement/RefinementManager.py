from pysgpp.extensions.datadriven.uq.operations import balance
from pysgpp import (DataVector, HashGridPoint,
                    SurplusRefinementFunctor,
                    HashGridStorage,
                    GridType_Linear, GridType_LinearL0Boundary,
                    GridType_LinearBoundary,
                    GridType_ModLinear,
                    GridType_Poly,
                    GridType_PolyBoundary)

import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import copyGrid, \
    polyGridTypes, bsplineGridTypes, linearGridTypes
from pysgpp.extensions.datadriven.uq.refinement.AdmissibleSet import AdmissibleSparseGridNodeSet
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder


class RefinementManager(object):

    def __init__(self, admissibleSet=None, criterion=None,
                 localRefinementStrategy=None, red=None,
                 maxLevel=30, verbose=False):
        """
        Constructor
        @param admissibleSet:
        @param criterion:
        @param localRefinementStrategy:
        @param red:
        @param maxLevel:
        @param verbose:
        """
        if not red:
            def red(v):
                return np.max(v, axis=1)
            self._red = red
        else:
            self._red = red

        self._adaptPoints = 0
        self._adaptRate = 0
        self._adaptThreshold = 0.0
        self._maxLevel = min(31, maxLevel)
        self._adaptTimeWindow = []
        self._balancing = False

        self._admissibleSet = admissibleSet
        self._criterion = criterion
        self._localRefinementStrategy = localRefinementStrategy
        self.refOnBorder = True
        self._averageWeightening = False

        self.verbose = verbose

    def setAdaptThreshold(self, value):
        self._adaptThreshold = value

    def setAdaptPoints(self, value):
        self._adaptPoints = value

    def setAdaptRate(self, value):
        self._adaptRate = value

    def setAdaptTimeWindow(self, window):
        self._adaptTimeWindow = window

    def getAdaptTimeWindow(self):
        return self._adaptTimeWindow

    def setAdaptMaxLevel(self, level):
        self._maxLevel = level

    def getAdaptMaxLevel(self):
        return self._maxLevel

    def hasBalancing(self):
        return self._balancing

    def setBalancing(self, balancing):
        self._balancing = balancing

    def hasAverageWeightening(self):
        return self._averageWeightening

    def setAverageWeightening(self, value):
        self._averageWeightening = value

    def setAdmissibleSetCreator(self, admissibleSet):
        self._admissibleSet = admissibleSet

    def getAdmissibleSet(self):
        return self._admissibleSet

    def setRefinementCriterion(self, criterion):
        self._criterion = criterion

    def getRefinementCriterion(self):
        return self._criterion

    def setLocalRefinementStrategy(self, localRefinementStrategy):
        self._localRefinementStrategy = localRefinementStrategy

    def getLocalRefinementStrategy(self):
        return self._localRefinementStrategy

    def refineOnTheBorder(self, refOnBorder):
        self.refOnBorder = refOnBorder

    # # Calculates the number of points which should be refined
    #
    # @param refinablePoints: integer number of points which can be refined
    # @return: integer number of point which should be refined
    def getNumOfPointsToRefine(self, refinablePoints):
        ratePoints = max(1, int(np.round(self._adaptRate * refinablePoints)))
        if self._adaptPoints == 0:
            return ratePoints
        elif self._adaptRate == 0:
            return self._adaptPoints
        else:
            return min(ratePoints, self._adaptPoints)


    def refineGrid(self, grid, knowledge, params=None, qoi="_", refinets=[0]):
        # check if this method is used in the right context
        if grid.getType() not in polyGridTypes + bsplineGridTypes + linearGridTypes:
                raise AttributeError('Grid type %s is not supported' %
                                     grid.getType())

        if params is None:
            # define standard uniform params
            uncertainParams = ParameterBuilder().defineUncertainParameters()
            for idim in xrange(grid.getStorage().getDimension()):
                uncertainParams.new().isCalled("x_%i" % idim).withUniformDistribution(0, 1)
            params = uncertainParams.andGetResult()

        # get refinement candidates
        if self.verbose:
            print "compute ranking"
        B = self.candidates(grid, knowledge, params, qoi, refinets)

        # now do the refinement
        return self.__refine(grid, B, simulate=False)


    def candidates(self, grid, knowledge, params, qoi="_", ts=None):
        """
        Load the candidates for refinement
        @param grid: Grid
        @param knowledge: ASGCKnowledge
        @param params: Parameter containing the marginal distributions
        @param qoi: string, quantity of interest
        @param ts: time steps of interest
        """
        # get knowledge type that is needed by refinement criterion
        dtype = self._criterion.getKnowledgeType()

        # get admissible set
        data = self._admissibleSet.values()

        v = np.ndarray([len(data), len(ts)], dtype='float')

        # run over the coefficients for all time steps
        for i, t in enumerate(ts):
            # get surpluses
            alphas = knowledge.getAlpha(qoi, t, dtype)
            # rank each admissible point
            for j, gp in enumerate(data):
                # run over all time steps for current grid point
                v[j, i] = self._criterion.rank(grid, gp, alphas, params, t)

#             # get the result
#             r = np.zeros(self._admissibleSet.getSize())
#             for i, gp in enumerate(self._admissibleSet.values()):
#                 r[i] = self._criterion.rank(grid, gp, alphas, params)
#             rix = np.argsort(r)
#             rjx = np.argsort(rix)
#             p = DataVector(2)
#             fig = plt.figure()
#             plotDensity2d(params.getIndependentJointDistribution())
#             values = self._admissibleSet.values()
#             for i, ri in enumerate(rix):
#                 gs.getCoordinates(values[ri], p)
#                 while i < len(rix) and.getStandardCoordinate(r[ri] - r[rix[i]]) < 1e-10:
#                     i += 1
#                 plt.plot(p[0], p[1], marker="o", color='yellow')
#                 plt.text(p[0], p[1], "%i" % (i - 1,), color='yellow',
#                          fontsize=12)
#             plt.xlim(0, 1)

#         # apply a subset ranking if the adding algorithm is used
#         if self._averageWeightening and \
#                 isinstance(self._admissibleSet, AdmissibleSparseGridNodeSet):
#             # simulate the refinement and sum over all subgroup values
#             # refine all the points in the admissible set
#             newGridPoints = {}
#             for j, gp in enumerate(data):
#                 newGridPoints[j] = self.__refine(learner, [(v[j, i], gp)],
#                                                  simulate=True)
#             for i, t in enumerate(ts):
#                 if self.verbose:
#                     print "compute merged ranking: %i/%i" % (i + 1, len(ts))
#
#                 # get surpluses
#                 alphas = learner.getKnowledge().getAlpha(qoi, t, dtype)
#                 # rank each admissible point
#                 s = 0.
#                 # take the refined points for gp in row j
#                 for j, points in newGridPoints.items():
#                     # sum up all the rankings
#                     for gp in points:
#                         key = (t, gp.getHash())
#                         if key not in d:
#                             d[key] = self._criterion.rank(grid, gp, alphas, params, t)
#                         s += d[key]
#                     # and take the mean as ranking
#                     v[j, i] = s / len(newGridPoints)

#             # simulate the refinement and sum over all subgroup values
#             for i, t in enumerate(ts):
#                 print "compute merged ranking: %i/%i" % (i + 1, len(ts))
#                 # get surpluses
#                 alphas = learner.getKnowledge().getAlpha(qoi, t, dtype)
#                 # rank each admissible point
#                 for j, gp in enumerate(data):
#                     newGridPoints = self.__refine(learner, [(v[j, i], gp)],
#                                                   simulate=True)
#                     s = 0.
#                     for gp in newGridPoints:
#                         s += self._criterion.rank(grid, gp, alphas, params)
#                     v[j, i] = s / len(newGridPoints)

        # reduce the collected contributions
        w = self._red(v)

        # get just the ones which are higher than the threshold
        B = [(y, data[ix]) for ix, y in enumerate(w) if y > self._adaptThreshold]

        # check if there are any nodes left
        if len(B) == 0:
            return []

        # sort by ranking
        B.sort(key=lambda x: x[0])

        return B

    def __refine(self, grid, B, simulate=False):
        # get sparse grid
        if simulate:
            grid = copyGrid(grid)

        # find how many points should be refined
        pointsNum = self.getNumOfPointsToRefine(len(B))

        # refine now step by step
        newGridPoints = []
        refinedPoints = []
        gs = grid.getStorage()

        # size of grid before refinement
        n1 = gs.getSize()

        # as long as the end of learning has not been reached, continue...
        while pointsNum > 0 and len(B) > 0:
            # note: the highest rated grid point is at the end of B
            vi, gp = B.pop()

            # some printing
            if not simulate and self.verbose:
                print "refine %i/%i (%i, %i) = %g" % \
                    (pointsNum, len(B) + 1, len(newGridPoints),
                     len(refinedPoints), vi)

            # refine the grid
            oldSize = grid.getSize()
            nps = self._localRefinementStrategy.refine(grid, gp)
            assert grid.getSize() == oldSize + len(nps)

            # ## set surplus vector such that just the desired point
            # ## is going to be refined and nothing else
            # oldgs = HashGridStorage(gs)
            # alpha = DataVector(gs.getSize())
            # alpha.setAll(0.0)
            # alpha[gs.getSequenceNumber(gp)] = 2.0
            # refFunc = SurplusRefinementFunctor(alpha, 1, 1)
            # ## TODO: try refineMaxLevel(refFunc, maxLevel)
            # grid.getGenerator().refine(refFunc)

            # nps = []
            # for i in xrange(gs.getSize()):
            #     if not oldgs.isContaining(gs.getPoint(i)):
            #         nps.append(i)

            # check there have been added some new points
            if len(nps) > 0:
                # if something has been refined then reduce the number
                # of points which should still be refined
                pointsNum -= 1

                # store which point has been refined
                refinedPoints.append(HashGridPoint(gp))
                newGridPoints += nps

        # balance the grid
        if self._balancing:
            newGridPoints += balance(grid)

#         if not simulate:
#             import matplotlib.pyplot as plt
#             gs = grid.getStorage()
#             p = DataVector(gs.getDimension())
#
#             fig = plt.figure()
#             for gp in refinedPoints:
#                 gs.getCoordinates(gp, p)
#                 plt.plot(p[0], p[1], marker='o', markersize=20,
#                          linestyle='', color='green')
#
#             for i in xrange(gs.getSize()):
#                 gpi = gs.getPoint(i)
#                 gs.getCoordinates(gpi, p)
#                 if gpi in self._admissibleSet:
#                     plt.plot(p[0], p[1], marker='o', markersize=10,
#                              linestyle='', color='orange')
#                 else:
#                     plt.plot(p[0], p[1], marker='o', markersize=10,
#                              linestyle='', color='blue')
#
#             for gp in newGridPoints:
#                 gs.getCoordinates(gp, p)
#                 plt.plot(p[0], p[1], marker='o', markersize=10,
#                          linestyle='', color='red')
#
#             plt.title("size = %i" % gs.getSize())
#             plt.xlim(0, 1)
#             plt.ylim(0, 1)
# #             plt.show()
#             plt.savefig('/home/franzefn/Desktop/tmp/var_atan_i%i.png' % gs.getSize())
#             plt.close(fig)

        # update admissible set
        if not simulate:
            self._admissibleSet.update(grid, newGridPoints)

        # make sure that I have collected all the new grid points
        assert len(newGridPoints) == gs.getSize() - n1

        return newGridPoints
