'''
Created on Feb 6, 2015

@author: franzefn
'''
from bin.uq.operations.sparse_grid import checkPositivity, \
    insertHierarchicalAncestors, insertPoint, getHierarchicalAncestors, copyGrid, \
    dehierarchize, hierarchize, hasChildren, countChildren
from pysgpp import GridIndex, createOperationEval, DataVector, IndexList
import warnings
from bin.uq.uq_plot.plot2d import plotSG2d
import pylab as plt
from interpolateParents import InterpolateParents


class OperationMakePositive(object):

    def __init__(self, grid):
        self.grid = grid

    def setInterpolationAlgorithm(self, algorithm):
        self.algorithm = algorithm

    def makeCurrentGridPositive(self, grid, alpha):
        nodalValues = dehierarchize(grid, alpha)
        cnt = 0
        for i, yi in enumerate(nodalValues.array()):
            if yi < 0:
                nodalValues[i] = 0
                cnt += 1
        if cnt > 0:
            alpha = hierarchize(grid, nodalValues)
            warnings.warn("negative function values encountered, this should not happen")

        return alpha

    def addFullGridPoints(self, grid, alpha):
        """
        Add all those full grid points with |l|_1 <= n, where n is the
        maximun level of the sparse grid
        @param grid: Grid sparse grid to be discretized
        @param alpha: DataVector hierarchical coefficients
        """
        # 1. create a new grid to work on
        gs = grid.getStorage()
        newGridPoints = []
        newGrid = copyGrid(grid)
        p = DataVector(gs.dim())
        opEval = createOperationEval(grid)

        # run over all grid points...
        for i in xrange(gs.size()):
            gp = gs.get(i)
            # ... and check if they are leaf nodes
            if countChildren(grid, gp) < 2 * gs.dim():
                # load all hierarchical ancestors
                ancestors = getHierarchicalAncestors(grid, gp)
                foundNegativeCoefficient = alpha[i] < 0
                j = 0
                # if there is at least one hierarchical ancestor that is
                # negative, then all full grid points in the support of the
                # current node can be negative...
                while not foundNegativeCoefficient and j < len(ancestors):
                    ix = gs.seq(ancestors[j][1])
                    foundNegativeCoefficient = alpha[ix] < 0
                    j += 1

                if foundNegativeCoefficient:
                    # ... therefore we run over all children up to the current
                    # max level and check if the function value is indeed
                    # negative
                    for d in xrange(gs.dim()):
                        if gp.getLevel(d) < gs.getMaxLevel():
                            children = [gp]
                            while len(children) > 0:
                                gpc = children.pop()
                                gpl = GridIndex(gpc)
                                gs.left_child(gpl, d)
                                gpr = GridIndex(gpc)
                                gs.right_child(gpr, d)

                                # if the function value for the left child
                                # is negtive, then add it with all its
                                # hierarchical ancestors
                                gpl.getCoords(p)
                                if opEval.eval(alpha, p) < 0:
                                    newGridPoints += insertPoint(newGrid, gpl)
                                    newGridPoints += insertHierarchicalAncestors(newGrid, gpl)
                                # if the function value for the right child
                                # is negtive, then add it with all its
                                # hierarchical ancestors
                                gpr.getCoords(p)
                                if opEval.eval(alpha, p) < 0:
                                    newGridPoints += insertPoint(newGrid, gpr)
                                    newGridPoints += insertHierarchicalAncestors(newGrid, gpr)

                                # store them for next round
                                if gpl.getLevel(d) < gs.getMaxLevel():
                                    children.append(gpl)
                                    children.append(gpr)

        # recompute the leaf property and return the result
        newGrid.getStorage().recalcLeafProperty()
        return newGrid, newGridPoints

    def coarsening(self, grid, alpha, newGridPoints):
        """
        Removes all unnecessary grid points. A grid point is defined as
        unnecessary if it is a leaf node and its hierarchical coefficient is
        negative. This is applied just for the newly added points.
        @param grid: Grid
        @param alpha: DataVector hierarchical coefficients
        @param newGridPoints: newly added grid points
        """
        newGrid = copyGrid(grid)
        newGs = newGrid.getStorage()
        newAlpha = DataVector(alpha)
        notAffectedGridPoints = []
        p = DataVector(newGs.dim())
        # remove all entirely negative weighted sub-branches in the tree
        # just consider the newly added ones by discretization
        while True:
            notAffectedGridPoints = []
            toBeRemoved = IndexList()
            for gp in newGridPoints:
                ix = newGs.seq(gp)
                gp.getCoords(p)
                # if the grid point is a leaf and has negative weight
                # we dont need it to make the function positive
                if hasChildren(grid, gp) and newAlpha[ix] < 0.:
                    toBeRemoved.append(ix)
                else:
                    notAffectedGridPoints.append(gp)

            # remove the identified grid points
            if toBeRemoved.size() > 0:
                newGs.deletePoints(toBeRemoved)
                newGs.recalcLeafProperty()
                # reset lists
                newGridPoints = notAffectedGridPoints
                # copy remaining alpha values
                newAlpha = DataVector(newGs.size())
                gs = grid.getStorage()
                for i in xrange(newGs.size()):
                    gp = newGs.get(i)
                    newAlpha[i] = alpha[gs.seq(gp)]
            else:
                break

        return newGrid, newAlpha, newGridPoints

    def makePositive(self, alpha):
        """
        insert recursively all grid points such that the function is positive
        defined. Interpolate the function values for the new grid points using
        the registered algorithm.
        @param alpha: DataVector hierarchical coefficients
        """
        # make sure that the function is positive at every existing grid point
        grid = self.grid
        alpha = self.makeCurrentGridPositive(grid, alpha)
        # start adding points
        newGridPoints = []
        while True:
            # add all those grid points which have an ancestor with negative
            # coefficient
            newGrid, addedGridPoints = self.addFullGridPoints(grid, alpha)

            if len(addedGridPoints) == 0:
                break

            newGridPoints += addedGridPoints
            nodalValues = dehierarchize(grid, alpha)
            nodalValues.resizeZero(newGrid.getSize())
            newAlpha = hierarchize(newGrid, nodalValues)
            newAlpha = self.algorithm.computeHierarchicalCoefficients(newGrid,
                                                                      newAlpha,
                                                                      addedGridPoints)

#             fig = plt.figure()
#             plotSG2d(newGrid, newAlpha)
#             plt.title("after density estimation")
#             fig.show()

            # check manually if all nodal points are positive
            newNodalValues = dehierarchize(newGrid, newAlpha)
            # collect all negative function values
            forceToBePositive = []
            newGs = newGrid.getStorage()
            for gp in addedGridPoints:
                if newNodalValues[newGs.seq(gp)] < 0.:
                    forceToBePositive.append(gp)

            if len(forceToBePositive) > 0:
                # if not, interpolate the function values for the negative
                # grid points
                warnings.warn("manually forcing the the function to be positive")
                newAlpha = self.makeCurrentGridPositive(newGrid, newAlpha)
                # newAlpha = InterpolateParents().computeHierarchicalCoefficients(newGrid, newAlpha, forceToBePositive)

#                 fig = plt.figure()
#                 plotSG2d(newGrid, newAlpha)
#                 plt.title("after setting negative grid points to zero")
#                 fig.show()

            grid, alpha = newGrid, newAlpha

        # coarsening: remove all new grid points with negative surplus
        coarsedGrid, coarsedAlpha, _ = self.coarsening(newGrid, newAlpha,
                                                       newGridPoints)

        # security check for positiveness
        checkPositivity(coarsedGrid, coarsedAlpha)

        return coarsedGrid, coarsedAlpha
