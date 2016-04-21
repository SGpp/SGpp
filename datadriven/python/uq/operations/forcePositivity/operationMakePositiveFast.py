'''
Created on Apr 19, 2016

@author: franzefn
'''
from pysgpp.extensions.datadriven.uq.operations import checkPositivity, \
    insertHierarchicalAncestors, insertPoint, copyGrid, \
    dehierarchize, hierarchize, hasChildren, hasAllChildren
from pysgpp import HashGridIndex, createOperationEval, DataVector, IndexList, \
    createOperationQuadrature, LinearBoundary, PolyBoundary
import warnings
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getHierarchicalAncestors, \
    insertTruncatedBorder, getBoundsOfSupport
import numpy as np
from matplotlib.patches import Rectangle
from pysgpp.extensions.datadriven.uq.transformation import LinearTransformation, \
    JointTransformation


class OperationMakePositiveFast(object):

    def __init__(self, grid, algorithm=None):
        self.grid = grid
        self.verbose = True
        self.algorithm = algorithm


    def plotDebug(self, newGrid, overlappingGridPoints):
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




    def setInterpolationAlgorithm(self, algorithm):
        self.algorithm = algorithm

    def makeCurrentNodalValuesPositive(self, grid, alpha):
        nodalValues = dehierarchize(grid, alpha)
        cnt = 0
        for i, yi in enumerate(nodalValues):
            if yi < 0:
                nodalValues[i] = 0
                cnt += 1
        if cnt > 0:
            alpha = hierarchize(grid, nodalValues)
#             if self.verbose:
#                 warnings.warn("negative function values at grid points encountered, this should not happen")

        return alpha

    
    def findGridPointsWithNegativeCoefficient(self, grid, alpha):
        negativeGridPoints = {}
        gs = grid.getStorage()
        for i in xrange(len(alpha)):
            if (alpha[i] < 0.0):
                negativeGridPoints[i] = gs.get(i)
        return negativeGridPoints


    def getOverlappingPoints(self, i, gpi, gps):
        numDims = gpi.getDimension()
        
        overlap = {}
        for j, gpj in gps.items():
            for d in xrange(numDims):
                # get level index
                lid, iid = gpi.getLevel(d), gpi.getIndex(d)
                ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

#                 # same level and different index
#                 if lid == ljd and iid != ijd:
#                     if (i, j) in overlap:
#                         del overlap[(i, j)]
#                     break

                # check if they have overlapping support
                xlowi, xhighi = getBoundsOfSupport(lid, iid)
                xlowj, xhighj = getBoundsOfSupport(ljd, ijd)
        
                xlow = max(xlowi, xlowj)
                xhigh = min(xhighi, xhighj)

                # different level but not ancestors
                if xlow >= xhigh:
                    if (i, j) in overlap:
                        del overlap[(i, j)]
                    break

                # they overlap in the current dimension
                if d == 0:
                    overlap[(i, j)] = (gpi, gpj)
        return overlap


    def findOverlappingBasisFunctions(self, grid, negativeGridPoints):
        overlappingGridPoints = {}
        
        gs = grid.getStorage()
        numDims = gs.getDimension()
        remainingGridPointsToCheck = negativeGridPoints.copy()
        for i, gpi in negativeGridPoints.items():
            # remove the current grid point from the list of comparison
            # due to symmetry
            del remainingGridPointsToCheck[i]

            overlap = self.getOverlappingPoints(i, gpi, remainingGridPointsToCheck)

            # copy the overlap to the return variable
            for key, value in overlap.items():
                overlappingGridPoints[key] = value

        return overlappingGridPoints


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


    def findIntersection(self, gpi, gpj):
        # find maximum level
        numDims = gpi.getDimension()
        gp = HashGridIndex(numDims)
        
        for d in xrange(numDims):
            if gpi.getLevel(d) > gpj.getLevel(d):
                gp.set(d, gpi.getLevel(d), gpi.getIndex(d))
            else:
                gp.set(d, gpj.getLevel(d), gpj.getIndex(d))
        return gp

    def findIntersections(self, overlappingGridPoints):
        intersections = [None] * len(overlappingGridPoints)
        for i, (gpi, gpj) in enumerate(overlappingGridPoints.values()):
            intersections[i] = self.findIntersection(gpi, gpj)
        return intersections
    
    
    def removeAlreadyExistingGridPoints(self, grid, intersections):
        gs = grid.getStorage()
        return [gp for gp in intersections if not gs.has_key(gp)]
    
    def addFullGridPoints(self, grid, alpha, intersections):
        """
        Add all those full grid points with |accLevel|_1 <= n, where n is the
        maximun level of the sparse grid
        @param grid: Grid sparse grid to be discretized
        @param intersections:
        """

        # make grid isotropic by adding all missing hierarchical ancestors
        newGridPoints = []
        # 3. insert them in a new grid if the function is negative
        p = DataVector(grid.getStorage().getDimension())
        opEval = createOperationEval(grid)
        alphaVec = DataVector(alpha)
        for gp in intersections:
            newGridPoints += insertPoint(grid, gp)
            newGridPoints += insertHierarchicalAncestors(grid, gp)

        # recompute the leaf property and return the result
        grid.getStorage().recalcLeafProperty()
        return newGridPoints


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

        newGridPoints = []
        while True:
            newAlpha = self.makeCurrentNodalValuesPositive(newGrid, newAlpha)

            if self.verbose:
                print "-" * 80
                print "# grid points      : %i" % (newGrid.getSize(),)

            negativeGridPoints = self.findGridPointsWithNegativeCoefficient(newGrid, newAlpha)

            if self.verbose:
                print "# negative coeffs  : %i/%i" % (len(negativeGridPoints), newGrid.getSize())

            overlappingGridPoints = self.findOverlappingBasisFunctions(newGrid, negativeGridPoints)

            if self.verbose:
                print "# overlapping basis: %i/%i" % (len(overlappingGridPoints), (len(negativeGridPoints) - 1) ** 2)

            intersections = self.findIntersections(overlappingGridPoints)

            if self.verbose:
                print "# intersections    : %i/%i" % (len(intersections), len(overlappingGridPoints))

            coarsedIntersections = self.removeAlreadyExistingGridPoints(newGrid, intersections)

            if self.verbose:
                print "# candidates       : %i/%i" % (len(coarsedIntersections), len(intersections))

            addedGridPoints = []
            if len(coarsedIntersections) > 0:
                addedGridPoints = self.addFullGridPoints(newGrid, newAlpha, coarsedIntersections)
                newAlpha = np.append(newAlpha, np.zeros(len(addedGridPoints)))
                newGridPoints += addedGridPoints

                if len(addedGridPoints) == 0:
                    break

                if self.verbose:
                    print "# new grid points  : %i" % len(addedGridPoints)

                if self.algorithm is not None:
                    # compute now the hierarchical coefficients for the newly added points
                    newAlpha = self.algorithm.computeHierarchicalCoefficients(newGrid,
                                                                              newAlpha,
                                                                              addedGridPoints)
            else:
                break

        # coarsening: remove all new grid points with zero surplus
        coarsedGrid, coarsedAlpha = self.coarsening(newGrid, newAlpha)
#         coarsedGrid, coarsedAlpha = newGrid, newAlpha
        if self.verbose:
            print "-" * 80
            print "# coarsed grid     : %i -> %i" % (newGrid.getSize(),
                                                       coarsedGrid.getSize(),)
            print "# full grid        : %i" % (2 ** self.grid.getStorage().getMaxLevel() - 1) ** self.grid.getStorage().getDimension()

        # security check for positiveness
        checkPositivity(coarsedGrid, coarsedAlpha)

        return coarsedGrid, coarsedAlpha


