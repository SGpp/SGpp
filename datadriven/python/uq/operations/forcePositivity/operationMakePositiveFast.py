'''
Created on Apr 19, 2016

@author: franzefn
'''
from pysgpp.extensions.datadriven.uq.operations import checkPositivity, \
    insertHierarchicalAncestors, insertPoint, copyGrid, \
    dehierarchize, hierarchize, hasChildren, hasAllChildren
from pysgpp import HashGridPoint, createOperationEval, DataVector, IndexList
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getHierarchicalAncestors, \
    insertTruncatedBorder, getBoundsOfSupport, evalSGFunctionMulti, \
    evalSGFunction
import numpy as np
from matplotlib.patches import Rectangle
from pysgpp.extensions.datadriven.uq.transformation import LinearTransformation, \
    JointTransformation
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.fullGridSearch import FullGridCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.localFullGridSearch import LocalFullGridCandidates
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.findIntersections import IntersectionCandidates


class OperationMakePositiveFast(object):

    def __init__(self, grid,
                 candidateSearchAlgorithm=None,
                 interpolationAlgorithm=None,):
        self.grid = grid
        self.numDims = grid.getStorage().getDimension()
        self.maxLevel = grid.getStorage().getMaxLevel()
        self.interpolationAlgorithm = interpolationAlgorithm
        self.candidateSearchAlgorithm = candidateSearchAlgorithm
        if self.candidateSearchAlgorithm is None:
            self.candidateSearchAlgorithm = IntersectionCandidates()

        self.maxNewGridPoints = 10
        self.addAllGridPointsOnNextLevel = True

        self.lastMinimumCandidateLevelSum = None

        self.verbose = True


    def plotDebugIntersections(self, newGrid, overlappingGridPoints):
        # -----------------------------------------------------------------
        # plot result
        gs = newGrid.getStorage()
        for n, ((i, j), (gpi, gpj)) in enumerate(overlappingGridPoints.items()):
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

            xlim = self.computeBoundsOfOverlappingPatch(gpi, gpj)
            currentAxis = plt.gca()
            currentAxis.add_patch(Rectangle((xlim[0, 0], xlim[1, 0]),
                                            np.diff(xlim[0, :2]), np.diff(xlim[1, :2]),
                                            facecolor="gray", alpha=0.9))

            gp = self.findIntersection(gpi, gpj)
            plt.plot(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1), "o ", color="black")

            plt.title("%i/%i" % (n + 1, len(overlappingGridPoints)))
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.show()
            plt.close(fig)


    def plotDebug(self, grid, alpha, addedGridPoints, candidates):
        # -----------------------------------------------------------------
        # plot result
        fig = plt.figure()
        plotSG2d(grid, alpha, show_negative=True, show_grid_points=True)
        plt.title("iteration = %i" % self.candidateSearchAlgorithm.iteration)
        
        p = DataVector(grid.getStorage().getDimension())
        for gp in candidates:
            gp.getStandardCoordinates(p)
            plt.plot(p[0], p[1], "o ", color="green")

        for gp in addedGridPoints:
            gp.getStandardCoordinates(p)
            plt.plot(p[0], p[1], "o ", color="yellow")

        fig.show()

    def setInterpolationAlgorithm(self, algorithm):
        self.interpolationAlgorithm = algorithm


    def setCandidateSetSearchAlgorithm(self, algorithm):
        self.candidateSearchAlgorithm = algorithm


    def makeAddedNodalValuesPositive(self, grid, alpha, addedGridPoints, tol=-1e-14):
        neg = []
        gs = grid.getStorage()
        x = DataVector(gs.getDimension())
        for gp in addedGridPoints:
            gp.getStandardCoordinates(x)
            yi = evalSGFunction(grid, alpha, x.array())
            if yi < tol:
                i = gs.getSequenceNumber(gp)
                alpha[i] -= yi
                assert alpha[i] > -1e-14
                assert evalSGFunction(grid, alpha, x.array()) < 1e-14
        return alpha


    def makeCurrentNodalValuesPositive(self, grid, alpha, tol=-1e-14):
        nodalValues = dehierarchize(grid, alpha)
        neg = []
        for i, yi in enumerate(nodalValues):
            if yi < tol:
                nodalValues[i] = 0
                neg.append(i)
        if len(neg) > 0:
            alpha = hierarchize(grid, nodalValues)

        return alpha
    

    def sortCandidatesByLevelSum(self, candidates):
        gps = {}
        for gp in candidates:
            levelSum = gp.getLevelSum()
            if levelSum not in gps:
                gps[levelSum] = [gp]
            else:
                gps[levelSum].append(gp)
        return gps


    def addFullGridPoints(self, grid, alpha, candidates, tol=-1e-14):
        """
        Add all those full grid points with |accLevel|_1 <= n, where n is the
        maximun level of the sparse grid
        @param grid: Grid sparse grid to be discretized
        @param candidates:
        @param tol:
        """
        # remove all the already existing candidates
        gs = grid.getStorage()
        nonExistingCandidates = [gp for gp in candidates if not gs.isContaining(gp)]

        # sort the non existing grid points by level
        finalCandidates = self.sortCandidatesByLevelSum(nonExistingCandidates)

        levelSums = sorted(finalCandidates.keys())
        if self.lastMinimumCandidateLevelSum is None:
            ix = 0
        else:
            ixs = np.where(np.array(levelSums) == self.lastMinimumCandidateLevelSum)[0]
            if len(ixs) == 0:
                ix = 0
            elif len(ixs) == 1:
                ix = ixs[0] + 1
            else:
                raise AttributeError("this should never happen")

        done = False
        addedGridPoints = []
        minLevelSum = -1
        nextLevelCosts = 0

        while not done and ix < len(levelSums):
            minLevelSum = levelSums[ix]
            self.lastMinimumCandidateLevelSum = minLevelSum

            currentCandidates = finalCandidates[minLevelSum]
            if self.verbose:
                print "# check candidates    : %i/%i (at |l|_1 = %i <= %i)" % (len(currentCandidates),
                                                                               len(candidates),
                                                                               minLevelSum,
                                                                               np.max(levelSums)),

            # evaluate the remaining candidates
            samples = np.ndarray((len(currentCandidates), self.numDims))
            p = DataVector(self.numDims)
            for j, gp in enumerate(currentCandidates):
                gp.getStandardCoordinates(p)
                samples[j, :] = p.array()
            eval = evalSGFunctionMulti(grid, alpha, samples)

            negativeNonExistingCandidates = [gp for j, gp in enumerate(currentCandidates) if eval[j] < tol]

            if self.verbose:
                print "-> %i : considered candidates" % len(negativeNonExistingCandidates)

            for gp in negativeNonExistingCandidates:
                addedGridPoints += insertPoint(grid, gp)
                addedGridPoints += insertHierarchicalAncestors(grid, gp)

                if not self.addAllGridPointsOnNextLevel and len(addedGridPoints) > self.maxNewGridPoints:
                    done = True
                    break

            if self.addAllGridPointsOnNextLevel and len(addedGridPoints) > 0:
                nextLevelCosts = len(finalCandidates[minLevelSum])
                done = True

            ix += 1

        # recompute the leaf property and return the result
        grid.getStorage().recalcLeafProperty()

        return addedGridPoints, minLevelSum, nextLevelCosts


    def coarsening(self, grid, alpha, newGridPoints, tol=1e-14):
        """
        Removes all unnecessary grid points. A grid point is defined as
        unnecessary if it is a leaf node and its hierarchical coefficient is
        negative. This is applied just for the newly added points.
        @param grid: Grid
        @param alpha: numpy array hierarchical coefficients
        """
        notAffectedGridPoints = []
        # remove all entirely negative weighted sub-branches in the tree
        # just consider the newly added ones by discretization
        iteration = 1
        while True:
            gs = grid.getStorage()

            notAffectedGridPoints = []
            toBeRemoved = IndexList()
            for gp in newGridPoints:
                # if the grid point is a leaf and has negative weight
                # we dont need it to make the function positive
                if gs.isContaining(gp):
                    ix = gs.getSequenceNumber(gp)
                    if gs.getPoint(ix).isLeaf() and np.abs(alpha[ix]) < tol:
                        toBeRemoved.append(ix)

            # remove the identified grid points
            if toBeRemoved.size() > 0:
                newGrid = copyGrid(grid)
                newGs = newGrid.getStorage()
                # delete the unneccessary grid points from the new grid
                newGs.deletePoints(toBeRemoved)
                newGs.recalcLeafProperty()
                # copy the remaining alpha values
                newAlpha = np.ndarray(newGs.getSize())
                for i in xrange(newGs.getSize()):
                    newAlpha[i] = alpha[gs.getSequenceNumber(newGs.getPoint(i))]

                grid, alpha = newGrid, newAlpha
                iteration += 1
            else:
                break

        return grid, alpha


    def makePositive(self, alpha):
        """
        insert recursively all grid points such that the function is positive
        defined. Interpolate the function values for the new grid points using
        the registered algorithm.
        @param alpha: numpy array hierarchical coefficients
        """
        # make sure that the function is positive at every existing grid point
        newGrid = copyGrid(self.grid)
        newAlpha = alpha

        if self.verbose:
            print

        iteration = 0
        newAlpha = self.makeCurrentNodalValuesPositive(newGrid, newAlpha)

        newGs = newGrid.getStorage()
        numDims = newGs.getDimension()
        newGridPoints = []
        addedGridPoints = {}

        numFullGridPoints = (2 ** self.maxLevel - 1) ** numDims
        minLevelSum = self.maxLevel + self.numDims - 1
        maxLevelSum = self.numDims * self.maxLevel
        totalCosts = 0
        currentCosts = 0
        candidates = []
        while minLevelSum < maxLevelSum and self.candidateSearchAlgorithm.hasMoreCandidates(newGrid, newAlpha, addedGridPoints):
            candidates, costs = self.candidateSearchAlgorithm.nextCandidateSet()
            totalCosts += costs
            oldGridSize = newGs.getSize()

            if self.verbose:
                print "iteration             : %i" % self.candidateSearchAlgorithm.iteration
                print "# found candidates    : %i/%i (costs = %i)" % (len(candidates), numFullGridPoints, costs)

            if len(candidates) > 0:
                addedGridPoints, minLevelSum, nextLevelCosts = self.addFullGridPoints(newGrid, newAlpha, candidates)
                currentCosts += nextLevelCosts
                assert oldGridSize + len(addedGridPoints) == newGs.getSize()
                if self.verbose:
                    print "# new grid points     : %i -> %i -> %i" % (oldGridSize,
                                                                      len(addedGridPoints),
                                                                      newGrid.getSize())
                    print "  current total costs : %i <= %i <= %i" % (currentCosts, totalCosts, numFullGridPoints)
                    print "-" * 80

                newAlpha = np.append(newAlpha, np.zeros(len(addedGridPoints)))

                if self.interpolationAlgorithm is not None:
                    # compute now the hierarchical coefficients for the newly added points
                    newAlpha = self.interpolationAlgorithm.computeHierarchicalCoefficients(newGrid,
                                                                                           newAlpha,
                                                                                           addedGridPoints)
                else:
                    newAlpha = self.makeAddedNodalValuesPositive(newGrid, newAlpha, addedGridPoints)

                # update list of new grid points
                newGridPoints += addedGridPoints

#                 if newGs.getDimension() == 2:
#                     self.plotDebug(newGrid, newAlpha, addedGridPoints, candidates)
            else:
                break

#         # coarsening: remove all new grid points with zero surplus
#         coarsedGrid, coarsedAlpha = self.coarsening(newGrid, newAlpha, newGridPoints)
#         if self.verbose:
#             print "                        old | coarsed | new | max (cand) | full"
#             print "# final grid          : %i <= %i <= %i <= %i (%i) <= %i" % (self.grid.getSize(),
#                                                                                coarsedGrid.getSize(),
#                                                                                newGrid.getSize(),
#                                                                                self.grid.getSize() + len(candidates),
#                                                                                len(candidates),
#                                                                                (2 ** self.maxLevel - 1) ** self.numDims)

        return newGrid, newAlpha
