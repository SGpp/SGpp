from pysgpp.extensions.datadriven.uq.operations import balance
from pysgpp import (DataVector, HashGridIndex,
                    SurplusRefinementFunctor,
                    HashGridStorage)

import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import copyGrid
from pysgpp.extensions.datadriven.uq.refinement.AdmissibleSet import AdmissibleSparseGridNodeSet


class Refinement(object):

    def __init__(self, admissibleSet=None, criterion=None,
                 localRefinementStrategy=None, red=None,
                 maxLevel=30):
        """
        Constructor
        @param admissibleSet:
        @param criterion:
        @param localRefinementStrategy:
        @param red:
        @param maxLevel:
        """
        if not red:
            def red(v):
                return np.max(v, axis=1)
            self._red = red
        else:
            self._red = red

        self._maxLevel = min(31, maxLevel)
        self._adaptTimeWindow = []
        self._balancing = False

        self._admissibleSet = admissibleSet
        self._criterion = criterion
        self._localRefinementStrategy = localRefinementStrategy
        self.refOnBorder = True
        self._averageWeightening = False

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

    def candidates(self, learner, ts=None):
        """
        Load the candidates for refinement
        @param learner: Learner
        @param ts: list of numeric, time steps
        """
        if ts is None or len(ts) == 0:
            # get time setting
            ts = learner.getTimeStepsOfInterest()

        qoi = learner.getQoI()
        grid = learner.getGrid()
        params = learner.getParameters()

        # get knowledge type that is needed by refinement criterion
        dtype = self._criterion.getKnowledgeType()

        # get admissible set
        data = self._admissibleSet.values()

        v = np.ndarray([len(data), len(ts)], dtype='float')

        d = {}

        # run over the coefficients for all time steps
        for i, t in enumerate(ts):
            # get surpluses
            alphas = learner.getKnowledge().getAlpha(qoi, t, dtype)
            # update refinement criterion
            self._criterion.update(grid, alphas, self._admissibleSet)
            # rank each admissible point
            for j, gp in enumerate(data):
                # run over all time steps for current grid point
                key = (t, gp.hash())
                if key not in d:
                    d[key] = self._criterion.rank(grid, gp, alphas, params)
                v[j, i] = d[key]

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
#                 values[ri].getCoords(p)
#                 while i < len(rix) and.getCoord(r[ri] - r[rix[i]]) < 1e-10:
#                     i += 1
#                 plt.plot(p[0], p[1], marker="o", color='yellow')
#                 plt.text(p[0], p[1], "%i" % (i - 1,), color='yellow',
#                          fontsize=12)
#             plt.xlim(0, 1)

        # apply a subset ranking if the adding algorithm is used
        if self._averageWeightening and \
                isinstance(self._admissibleSet, AdmissibleSparseGridNodeSet):
            # simulate the refinement and sum over all subgroup values
            # refine all the points in the admissible set
            newGridPoints = {}
            for j, gp in enumerate(data):
                newGridPoints[j] = self.__refine(learner, [(v[j, i], gp)],
                                                 simulate=True)
            for i, t in enumerate(ts):
                print "compute merged ranking: %i/%i" % (i + 1, len(ts))
                # get surpluses
                alphas = learner.getKnowledge().getAlpha(qoi, t, dtype)
                # rank each admissible point
                s = 0.
                # take the refined points for gp in row j
                for j, points in newGridPoints.items():
                    # sum up all the rankings
                    for gp in points:
                        key = (t, gp.hash())
                        if key not in d:
                            d[key] = self._criterion.rank(grid, gp, alphas, params)
                        s += d[key]
                    # and take the mean as ranking
                    v[j, i] = s / len(newGridPoints)

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
        x = learner.getAdaptThreshold()
        B = [(y, data[ix]) for ix, y in enumerate(w) if y > x]

        # check if there are any nodes left
        if len(B) == 0:
            return []

        # sort by ranking
        B.sort(key=lambda x: x[0])

        return B

    def refineGrid(self, learner, ts=None):
        # check if this method is used in the right context
        if learner.getGrid().getType() not in ('linear',
                                               'linearBoundary',
                                               'linearTruncatedBoundary',
                                               'modlinear',
                                               'ultraPoly',
                                               'ultraPolyTruncatedBoundary',
                                               'myPoly'):
            raise AttributeError('Grid type %s is not supported' %
                                 learner.getGrid().getType())

        # get refinement candidates
        print "compute ranking"
        B = self.candidates(learner, ts)

        # now do the refinement
        return self.__refine(learner, B, simulate=False)

    def __refine(self, learner, B, simulate=False):
        # get sparse grid
        grid = learner.getGrid()
        if simulate:
            oldGrid = grid
            grid = copyGrid(grid)
            learner.grid = grid

        # find how many points should be refined
        pointsNum = learner.getNumOfPointsToRefine(len(B))

        # refine now step by step
        newGridPoints = []
        refinedPoints = []
        gs = grid.getStorage()
        iteration = learner.iteration
        # size of grid before refinement
        n1 = gs.size()

        # as long as the end of learning has not been reached, continue...
        while pointsNum > 0 and len(B) > 0 and \
            (not learner.stopPolicy or
             not learner.stopPolicy.hasLimitReached(learner)):
            # note: the highest rated grid point is at the end of B
            vi, gp = B.pop()

            # some printing
            if not simulate:
                print "refine %i/%i (%i, %i) = %g" % \
                    (pointsNum, len(B), len(newGridPoints),
                     len(refinedPoints), vi)

            # refine the grid
            nps = self._localRefinementStrategy.refine(grid, gp)

            # ## set surplus vector such that just the desired point
            # ## is going to be refined and nothing else
            # oldgs = HashGridStorage(gs)
            # alpha = DataVector(gs.size())
            # alpha.setAll(0.0)
            # alpha[gs.seq(gp)] = 2.0
            # refFunc = SurplusRefinementFunctor(alpha, 1, 1)
            # ## TODO: try refineMaxLevel(refFunc, maxLevel)
            # grid.createGridGenerator().refine(refFunc)

            # nps = []
            # for i in xrange(gs.size()):
            #     if not oldgs.has_key(gs.get(i)):
            #         nps.append(i)

            # check there have been added some new points
            if not learner.stopPolicy or \
                    learner.stopPolicy.hasGridSizeChanged(learner):
                # if something has been refined then reduce the number
                # of points which should still be refined
                pointsNum -= 1

                # store which point has been refined
                refinedPoints.append(HashGridIndex(gp))
                newGridPoints += nps

                # increase iteration of the learner
                learner.iteration += 1

        # balance the grid
        if self._balancing:
            newGridPoints += balance(grid)

        # update admissible set
        if not simulate:
            self._admissibleSet.update(grid, newGridPoints)

        # make sure that I have collected all the new grid points
        assert len(newGridPoints) == gs.size() - n1

        # reset the iteration variable @TODO: the iteration variable
        # is ambiguous. It represents in the TrainingStopPolicy the
        # number of refinement steps, in the context of ASGC it
        # represents the number of refinements. So here we neglect the
        # first part and use it just internally so that the
        # hasLimitReached works.
        learner.iteration = iteration

#         if not simulate:
#             gs = grid.getStorage()
#             p = DataVector(gs.dim())
#
#             for gp in refinedPoints:
#                 gp.getCoords(p)
#                 plt.plot(p[0], p[1], marker='o', markersize=20,
#                          linestyle='', color='green')
#
#             for i in xrange(gs.size()):
#                 gs.get(i).getCoords(p)
#                 plt.plot(p[0], p[1], marker='o', markersize=10,
#                          linestyle='', color='blue')
#
#             for gp in newGridPoints:
#                 gp.getCoords(p)
#                 plt.plot(p[0], p[1], marker='o', markersize=10,
#                          linestyle='', color='red')
#
#             plt.title("size = %i" % gs.size())
#             plt.xlim(0, 1)
#             plt.ylim(0, 1)
#             plt.savefig('%i.png' % learner.iteration)

        # reset the learner if the refinement is just simulated
        if simulate:
            learner.grid = oldGrid

        return newGridPoints
