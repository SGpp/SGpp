'''
Created on Apr 19, 2016

@author: franzefn
'''
from pysgpp.extensions.datadriven.uq.operations import checkPositivity, \
    insertHierarchicalAncestors, insertPoint, copyGrid, \
    dehierarchize, hierarchize, hasChildren, hasAllChildren
from pysgpp import HashGridIndex, createOperationEval, DataVector, IndexList
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getHierarchicalAncestors, \
    insertTruncatedBorder, getBoundsOfSupport, evalSGFunctionMulti
import numpy as np
from matplotlib.patches import Rectangle
from pysgpp.extensions.datadriven.uq.transformation import LinearTransformation, \
    JointTransformation
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.fullGridSearch import FullGridCandidates


class OperationMakePositiveFast(object):

    def __init__(self, grid,
                 candidateSetAlgorithm=None,
                 interpolationAlgorithm=None,
                 candidateSearchAlgorithm=None):
        self.grid = grid
        self.numDims = grid.getStorage().getDimension()
        self.maxLevel = grid.getStorage().getMaxLevel()
        self.verbose = True
        self.interpolationAlgorithm = interpolationAlgorithm
        self.candidateSearchAlgorithm = candidateSearchAlgorithm
        if self.candidateSearchAlgorithm is None:
            self.candidateSearchAlgorithm = FullGridCandidates(grid)

        self.maxNewGridPoints = 10
        self.addAllGridPointsOnNextLevel = True


    def plotDebugIntersections(self, newGrid, overlappingGridPoints):
        # -----------------------------------------------------------------
        # plot result
        gs = newGrid.getStorage()
        for n, ((i, j), (gpi, gpj)) in enumerate(overlappingGridPoints.items()):
            fig = plt.figure()
            for k in xrange(gs.getSize()):
                gp = gs.get(k)
                x, y = gp.getCoord(0), gp.getCoord(1)
                if alpha[k] < 0.0:
                    plt.plot(x, y, "v ", color="red")
                else:
                    plt.plot(x, y, "^ ", color="white")

            # annotate the
            plt.annotate(str(i), (gpi.getCoord(0), gpi.getCoord(1)))
            plt.annotate(str(j), (gpj.getCoord(0), gpj.getCoord(1)))

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
            plt.plot(gp.getCoord(0), gp.getCoord(1), "o ", color="black")

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
            gp.getCoords(p)
            plt.plot(p[0], p[1], "o ", color="green")

        for gp in addedGridPoints:
            gp.getCoords(p)
            plt.plot(p[0], p[1], "o ", color="yellow")

        fig.show()

    def setInterpolationAlgorithm(self, algorithm):
        self.interpolationAlgorithm = algorithm


    def setCandidateSetSearchAlgorithm(self, algorithm):
        self.candidateSearchAlgorithm = algorithm


    def makeCurrentNodalValuesPositive(self, grid, alpha, addedGridPoints=None):
        nodalValues = dehierarchize(grid, alpha)
        neg = []
        for i, yi in enumerate(nodalValues):
            if yi < 0:
                nodalValues[i] = 0
                neg.append(i)
        if len(neg) > 0:
            alpha = hierarchize(grid, nodalValues)
#             if self.verbose:
#                 warnings.warn("negative function values at grid points encountered, this should not happen")

            # check if the coefficients of the new grid points are positive
            if addedGridPoints is not None:
                gs = grid.getStorage()
                assert all([alpha[gs.seq(gp)] >= -1e-13 for gp in addedGridPoints])
        return alpha

    

    def computeBoundsOfOverlappingPatch(self, gpi, gpj):
        # compute bounds of the overlapping patch
        numDims = gpi.getDimension()
        xlim = np.ndarray((numDims, 2))
        for d in xrange(numDims):
            # get level index
            lid, iid = gpi.getLevel(d), gpi.getIndex(d)
            ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

            # check if they have overlapping support
            xlowi, xhighi = getBoundsOfSupport(lid, iid)
            xlowj, xhighj = getBoundsOfSupport(ljd, ijd)

            xlim[d, 0] = max(xlowi, xlowj)
            xlim[d, 1] = min(xhighi, xhighj)

        return xlim


    def sortCandidatesByLevelSum(self, candidates):
        gps = {}
        for gp in candidates:
            levelSum = gp.getLevelSum()
            if levelSum not in gps:
                gps[levelSum] = [gp]
            else:
                gps[levelSum].append(gp)
        return gps


    def addFullGridPoints(self, grid, alpha, candidates):
        """
        Add all those full grid points with |accLevel|_1 <= n, where n is the
        maximun level of the sparse grid
        @param grid: Grid sparse grid to be discretized
        @param candidates:
        """
        # remove all the already existing candidates
        gs = grid.getStorage()
        nonExistingCandidates = [gp for gp in candidates if not gs.has_key(gp)]

        # evaluate the remaining candidates
        samples = np.ndarray((len(nonExistingCandidates), self.numDims))
        p = DataVector(self.numDims)
        for i, gp in enumerate(nonExistingCandidates):
            gp.getCoords(p)
            samples[i, :] = p.array()
        eval = evalSGFunctionMulti(grid, alpha, samples)

        negativeNonExistingCandidates = [gp for i, gp in enumerate(nonExistingCandidates) if eval[i] < 0.0]

        # sort the non existing grid points by level
        finalCandidates = self.sortCandidatesByLevelSum(negativeNonExistingCandidates)

        levelSums = sorted(finalCandidates.keys())
        done = False
        i = 0
        addedGridPoints = []
        minLevelSum = -1

#         for levelSum, gridPoints in gps.items():
#             print "# candidates at |l_i|_1 = %i: %i/%i" % (levelSum,
#                                                            len(gps[levelSum]),
#                                                            cnt)

        while not done and i < len(levelSums):
            minLevelSum = levelSums[i]
            for gp in finalCandidates[minLevelSum]:
                addedGridPoints += insertPoint(grid, gp)
                addedGridPoints += insertHierarchicalAncestors(grid, gp)

                if not self.addAllGridPointsOnNextLevel and len(addedGridPoints) > self.maxNewGridPoints:
                    done = True
                    break

            if self.addAllGridPointsOnNextLevel and len(addedGridPoints) > 0:
                done = True

            i += 1

        # recompute the leaf property and return the result
        grid.getStorage().recalcLeafProperty()

        return addedGridPoints, minLevelSum


    def coarsening(self, grid, alpha):
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
            for ix in xrange(gs.getSize()):
                # if the grid point is a leaf and has negative weight
                # we dont need it to make the function positive
                if gs.get(ix).isLeaf() and np.abs(alpha[ix]) < 1e-14:
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
                    newAlpha[i] = alpha[gs.seq(newGs.get(i))]

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
        addedGridPoints = {}

        numFullGridPoints = (2 ** self.maxLevel - 1) ** numDims
        minLevelSum = self.maxLevel + self.numDims - 1
        maxLevelSum = self.numDims * self.maxLevel
        totalCosts = 0
        while minLevelSum < maxLevelSum and self.candidateSearchAlgorithm.hasMoreCandidates(newGrid, newAlpha, addedGridPoints):
            candidates, costs = self.candidateSearchAlgorithm.nextCandidateSet()
            totalCosts += costs
            oldGridSize = newGs.getSize()

            if self.verbose:
                print "iteration             : %i" % self.candidateSearchAlgorithm.iteration
                print "# found candidates    : %i/%i (costs = %i)" % (len(candidates), numFullGridPoints, costs)

            addedGridPoints = {}
            if len(candidates) > 0:
                addedGridPoints, minLevelSum = self.addFullGridPoints(newGrid, newAlpha, candidates)
                assert oldGridSize + len(addedGridPoints) == newGs.getSize()
                if self.verbose:
                    print "# new grid points     : %i -> %i -> %i at |l|_1 = %i <= %i" % (oldGridSize,
                                                                                          len(addedGridPoints),
                                                                                          newGrid.getSize(),
                                                                                          minLevelSum,
                                                                                          maxLevelSum)
                    print "  current total costs : %i <= %i" % (totalCosts, numFullGridPoints)
                    print "-" * 80

                newAlpha = np.append(newAlpha, np.zeros(len(addedGridPoints)))

                if self.interpolationAlgorithm is not None:
                    # compute now the hierarchical coefficients for the newly added points
                    newAlpha = self.interpolationAlgorithm.computeHierarchicalCoefficients(newGrid,
                                                                                           newAlpha,
                                                                                           addedGridPoints)
                else:
                    newAlpha = self.makeCurrentNodalValuesPositive(newGrid, newAlpha, addedGridPoints)

#                 if newGs.getDimension() == 2:
#                     self.plotDebug(newGrid, newAlpha, addedGridPoints, candidates)
            else:
                break

        # coarsening: remove all new grid points with zero surplus
        coarsedGrid, coarsedAlpha = self.coarsening(newGrid, newAlpha)
        if self.verbose:
            print "# coarsed grid        : %i -> %i" % (newGrid.getSize(),
                                                        coarsedGrid.getSize())
            print "# full grid           :       %i" % (2 ** self.maxLevel - 1) ** self.numDims

        # security check for positiveness
        neg = checkPositivity(coarsedGrid, coarsedAlpha)
        
#         if len(neg) > 0:
#             # check at which grid points the function is negative
#             for i, (yi, gp) in neg.items():
#                     print "|%s|_1 = %i, %s -> %g" % ([gp.getLevel(d) for d in xrange(numDims)],
#                                                      np.sum([gp.getLevel(d) for d in xrange(numDims)]),
#                                                      [gp.getIndex(d) for d in xrange(numDims)],
#                                                      yi)

        return coarsedGrid, coarsedAlpha


