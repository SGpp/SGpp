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

    
    def findGridPointsWithNegativeCoefficient(self, gps, alpha):
        negativeGridPoints = {}
        for i, gp in gps.items():
            if alpha[i] < 0.0:
                negativeGridPoints[i] = gp
        return negativeGridPoints


    def separateGridPoints(self, gps, alpha):
        negativeGridPoints = {}
        positiveGridPoints = {}
        for i, gp in gps.items():
            if alpha[i] < 0.0:
                negativeGridPoints[i] = gp
            else:
                positiveGridPoints[i] = gp

        return negativeGridPoints, positiveGridPoints


    def findIntersectionsOfOverlappingSuppportsForOneGridPoint(self, i, gpi, gpsj, grid):
        numDims = gpi.getDimension()
        gs = grid.getStorage()
        overlap = {}
        already_checked = {}
        # find all possible intersections of grid points
        for j, gpj in gpsj.items():
            if i != j:
                level, index = self.findIntersection(gpi, gpj)
                key = tuple(level + index)
                if key not in overlap and key not in already_checked:
                    intersection = HashGridIndex(numDims)
                    for idim in xrange(numDims):
                        intersection.set(idim, level[idim], index[idim])

                    # check if the current grid point already exists
                    if not gs.has_key(intersection):
                        idim = 0
                        while idim < numDims:
                            # get level index
                            lid, iid = gpi.getLevel(idim), gpi.getIndex(idim)
                            ljd, ijd = gpj.getLevel(idim), gpj.getIndex(idim)

                            # same level and different index
                            if lid == ljd and iid != ijd:
                                break

                            # check if they have overlapping support
                            xlowi, xhighi = getBoundsOfSupport(lid, iid)
                            xlowj, xhighj = getBoundsOfSupport(ljd, ijd)

                            xlow = max(xlowi, xlowj)
                            xhigh = min(xhighi, xhighj)

                            # different level but not ancestors
                            if xlow >= xhigh:
                                break

                            idim += 1

                        # check whether the supports are overlapping
                        # in all dimensions
                        if idim == numDims:
                            overlap[key] = intersection
                else:
                    already_checked[key] = False
                
        return overlap


    def findIntersectionsOfOverlappingSuppports(self, gpsi, gpsj, grid):
        overlappingGridPoints = {}
        
        for i, gpi in gpsi.items():
            overlap = self.findIntersectionsOfOverlappingSuppportsForOneGridPoint(i, gpi, gpsj, grid)
            overlappingGridPoints.update(overlap)

        return overlappingGridPoints.values()


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
        level = [0] * numDims
        index = [0] * numDims
        
        for d in xrange(numDims):
            if gpi.getLevel(d) > gpj.getLevel(d):
                level[d] = gpi.getLevel(d)
                index[d] = gpi.getIndex(d)
            else:
                level[d] = gpj.getLevel(d)
                index[d] = gpj.getIndex(d)

        return level, index

    
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
        addedGridPoints = []
        # 3. insert them in a new grid if the function is negative
        p = DataVector(grid.getStorage().getDimension())
        opEval = createOperationEval(grid)
        alphaVec = DataVector(alpha)
        for gp in intersections:
            gp.getCoords(p)
            if opEval.eval(alphaVec, p) < 0.0:
                addedGridPoints += insertPoint(grid, gp)
                addedGridPoints += insertHierarchicalAncestors(grid, gp)
        # recompute the leaf property and return the result
        grid.getStorage().recalcLeafProperty()

        newGridPoints = {}
        gs = grid.getStorage()
        for gp in addedGridPoints:
            newGridPoints[gs.seq(gp)] = gp

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

        iteration = 0
        newAlpha = self.makeCurrentNodalValuesPositive(newGrid, newAlpha)

        newGs = newGrid.getStorage()
        addedGridPoints = dict([(i, newGs.get(i)) for i in xrange(newGs.getSize())])
        N0 = self.findGridPointsWithNegativeCoefficient(addedGridPoints, newAlpha)
        N1 = {}

        while True:
            cntChecks = len(N0) * len(addedGridPoints)
            if self.verbose:
                print "-" * 80
                print "# search intersections: %i = %i x %i" % (len(N0) * len(addedGridPoints),
                                                                len(N0), len(addedGridPoints))

            intersections = self.findIntersectionsOfOverlappingSuppports(N0, addedGridPoints, newGrid)

            if self.verbose:
                print "# found intersections : %i/%i" % (len(intersections),
                                                         cntChecks)

            addedGridPoints = []
            if len(intersections) > 0:
                addedGridPoints = self.addFullGridPoints(newGrid, newAlpha, intersections)

                if self.verbose:
                    print "# new grid points     : %i -> %i" % (len(addedGridPoints), newGrid.getSize())

                if len(addedGridPoints) == 0:
                    break

                newAlpha = np.append(newAlpha, np.zeros(len(addedGridPoints)))

                if self.algorithm is not None:
                    # compute now the hierarchical coefficients for the newly added points
                    newAlpha = self.algorithm.computeHierarchicalCoefficients(newGrid,
                                                                              newAlpha,
                                                                              addedGridPoints.values())
                else:
                    newAlpha = self.makeCurrentNodalValuesPositive(newGrid, newAlpha)
            else:
                break

            N0.update(N1)  # all the negative ones from last iteration
            N1 = self.findGridPointsWithNegativeCoefficient(addedGridPoints, newAlpha)

            iteration += 1

        # coarsening: remove all new grid points with zero surplus
        coarsedGrid, coarsedAlpha = self.coarsening(newGrid, newAlpha)
#         coarsedGrid, coarsedAlpha = newGrid, newAlpha
        if self.verbose:
            print "-" * 80
            print "# coarsed grid        : %i -> %i" % (newGrid.getSize(),
                                                     coarsedGrid.getSize())
            print "# full grid           :       %i" % (2 ** self.grid.getStorage().getMaxLevel() - 1) ** self.grid.getStorage().getDimension()

        # security check for positiveness
        checkPositivity(coarsedGrid, coarsedAlpha)

        return coarsedGrid, coarsedAlpha


