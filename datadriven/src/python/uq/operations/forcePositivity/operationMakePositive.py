'''
Created on Feb 6, 2015

@author: franzefn
'''
from bin.uq.operations import checkPositivity, \
    insertHierarchicalAncestors, insertPoint, copyGrid, \
    dehierarchize, hierarchize, hasChildren, hasAllChildren
from pysgpp import GridIndex, createOperationEval, DataVector, IndexList, \
    createOperationQuadrature
import warnings
from bin.uq.uq_plot.plot2d import plotSG2d
import pylab as plt
from bin.uq.operations.sparse_grid import getHierarchicalAncestors, \
    insertTrapezoidBorder


class OperationMakePositive(object):

    def __init__(self, grid):
        self.grid = grid
        self.verbose = False

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
            if self.verbose:
                warnings.warn("negative function values encountered, this should not happen")

        return alpha

    def findCandidatesSweep1d(self, d, gp, alpha, grid, acc,
                              negativeAncestorFound=False):
        gs = grid.getStorage()
        ix = gs.seq(gp)
        negativeAncestorFound |= alpha[ix] < 0.

        # just leaf nodes are possible candidates
        if hasAllChildren(grid, gp) and negativeAncestorFound:
            acc.append(gp)

        # get left child
        level, index = gp.getLevel(d), gp.getIndex(d)
        gs.left_child(gp, d)
        if gs.has_key(gp):
            self.findCandidatesSweep1d(d, gp, alpha, grid, acc,
                                       negativeAncestorFound)

        # get right child
        gp.set(d, level, index)
        gs.right_child(gp, d)
        if gs.has_key(gp):
            self.findCandidatesSweep1d(d, gp, alpha, grid, acc,
                                       negativeAncestorFound)
        gp.set(d, level, index)

    def findCandidates(self, grid, alpha):
        gs = grid.getStorage()
        candidates = {}

        # lookup dimension-wise
        for d in xrange(gs.dim()):
            # compute starting points by level sum
            anchors = []
            for i in xrange(gs.size()):
                accLevel = gs.get(i).getLevel(d)
                if accLevel == 1:
                    anchors.append(i)

            while len(anchors) > 0:
                # get next starting node
                ix = anchors.pop(0)
                gp = gs.get(ix)
                acc = []
                self.findCandidatesSweep1d(d, gp, alpha, grid, acc, False)
                # store candidates
                for gp in acc:
                    ix = gs.seq(gp)
                    if ix not in candidates:
                        candidates[ix] = GridIndex(gp)

        return candidates.values()

#     def findCandidates(self, grid, alpha):
#         # 1. create a new grid to work on
#         gs = grid.getStorage()
#         candidates = {}
#
#         # run over all grid points...
#         for i in xrange(gs.size()):
#             gp = gs.get(i)
#             # ... and check if they are leaf nodes
#             if not hasAllChildren(grid, gp):
#                 # load all hierarchical ancestors
#                 ancestors = getHierarchicalAncestors(grid, gp)
#                 foundNegativeCoefficient = alpha[i] < 0
#                 j = 0
#                 # if there is at least one hierarchical ancestor that is
#                 # negative, then all full grid points in the support of the
#                 # current node can be negative...
#                 while not foundNegativeCoefficient and j < len(ancestors):
#                     ix = gs.seq(ancestors[j][1])
#                     foundNegativeCoefficient = alpha[ix] < 0
#                     j += 1
#
#                 if foundNegativeCoefficient and i not in candidates:
#                     candidates[i] = gp
#
#         return candidates.values()

    def lookupFullGridPointsRec1d(self, grid, alpha, gp, d, p, opEval, maxLevel, acc):
        gs = grid.getStorage()
        level, index = gp.getLevel(d), gp.getIndex(d)

        # if the function value for the left child
        # is negtive, then add it with all its
        # hierarchical ancestors
        gs.left_child(gp, d)
        gp.getCoords(p)
        if opEval.eval(alpha, p) < 0:
            acc.append(GridIndex(gp))

        if level + 1 < maxLevel:
            self.lookupFullGridPointsRec1d(grid, alpha, gp, d, p, opEval, maxLevel, acc)

        # if the function value for the right child
        # is negtive, then add it with all its
        # hierarchical ancestors
        gp.set(d, level, index)
        gs.right_child(gp, d)
        gp.getCoords(p)
        if opEval.eval(alpha, p) < 0:
            acc.append(GridIndex(gp))

        # store them for next round
        if level + 1 < maxLevel:
            self.lookupFullGridPointsRec1d(grid, alpha, gp, d, p, opEval, maxLevel, acc)

        # reset the grid point
        gp.set(d, level, index)

    def lookupFullGridPoints(self, grid, alpha, candidates):
        acc = []
        gs = grid.getStorage()
        p = DataVector(gs.dim())
        opEval = createOperationEval(grid)
        # TODO: find local max level for adaptively refined grids
        maxLevel = gs.getMaxLevel()
        if "Boundary" in grid.getType():
            maxLevel += 1

        for gp in candidates:
            for d in xrange(gs.dim()):
                if 0 < gp.getLevel(d) < maxLevel:
                    self.lookupFullGridPointsRec1d(grid, alpha, gp, d, p,
                                                   opEval, maxLevel, acc)
        return acc

    def addFullGridPoints(self, grid, alpha):
        """
        Add all those full grid points with |accLevel|_1 <= n, where n is the
        maximun level of the sparse grid
        @param grid: Grid sparse grid to be discretized
        @param alpha: DataVector hierarchical coefficients
        """

        # make grid isotropic by adding all missing hierarchical ancestors
        newGridPoints = []
        # 1. locate candidates to be refined
        candidates = self.findCandidates(grid, alpha)
        # 2. do full grid search locally
        gridPoinsToBeAdded = self.lookupFullGridPoints(grid, alpha, candidates)
        # 3. insert them in a new grid
        for gp in gridPoinsToBeAdded:
            newGridPoints += insertPoint(grid, gp)
            newGridPoints += insertHierarchicalAncestors(grid, gp)
            if "Boundary" in grid.getType():
                newGridPoints += insertTrapezoidBorder(grid, gp)

        # recompute the leaf property and return the result
        grid.getStorage().recalcLeafProperty()
        return newGridPoints

    def coarsening(self, grid, alpha, newGridPoints):
        """
        Removes all unnecessary grid points. A grid point is defined as
        unnecessary if it is a leaf node and its hierarchical coefficient is
        negative. This is applied just for the newly added points.
        @param grid: Grid
        @param alpha: DataVector hierarchical coefficients
        @param newGridPoints: newly added grid points
        """
        gs = grid.getStorage()
        newAlpha = DataVector(alpha)
        notAffectedGridPoints = []
        p = DataVector(gs.dim())
        # remove all entirely negative weighted sub-branches in the tree
        # just consider the newly added ones by discretization
        iteration = 1
        while True:
            newGrid = copyGrid(grid)
            notAffectedGridPoints = []
            toBeRemoved = IndexList()
            for gp in newGridPoints:
                ix = gs.seq(gp)
                gp.getCoords(p)
                # if the grid point is a leaf and has negative weight
                # we dont need it to make the function positive
                if not hasChildren(grid, gp) and newAlpha[ix] < 0.:
                    toBeRemoved.append(ix)
                else:
                    notAffectedGridPoints.append(gp)

            # remove the identified grid points
            if toBeRemoved.size() > 0:
                # copy the old grid
                newGs = newGrid.getStorage()
                # delete the unneccessary grid points from the new grid
                newGs.deletePoints(toBeRemoved)
                newGs.recalcLeafProperty()
                # reset lists
                newGridPoints = notAffectedGridPoints
                # copy the remaining alpha values
                newAlpha = DataVector(newGs.size())
                for i in xrange(newGs.size()):
                    newAlpha[i] = alpha[gs.seq(newGs.get(i))]

                grid, gs, alpha = newGrid, newGs, newAlpha
                iteration += 1
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
        iteration = 1
        while True:
            # copy the old grid
            newGrid = copyGrid(grid)
            # add all those grid points which have an ancestor with negative
            # coefficient
            if self.verbose:
                print "-" * 60
                print "%i:" % iteration
                print "adding full grid points"
            addedGridPoints = self.addFullGridPoints(newGrid, alpha)

            if len(addedGridPoints) == 0:
                newAlpha = alpha
                break

            if self.verbose:
                print "learning the new density"

            newGridPoints += addedGridPoints
            # set the function value at the new grid points to zero
            newNodalValues = dehierarchize(grid, alpha)
            newNodalValues.resizeZero(newGrid.getSize())
            newAlpha = DataVector(alpha)
            newAlpha.resizeZero(newGrid.getSize())
            # compute now the hierarchical coefficients for the newly
            # added points
            newAlpha = self.algorithm.computeHierarchicalCoefficients(newGrid,
                                                                      newAlpha,
                                                                      addedGridPoints)
            # the function does not have to be positive now -> check it again
            if self.verbose:
                print "force function to be positive"
            # check manually if all nodal points are positive
            newNodalValues = dehierarchize(newGrid, newAlpha)
            # collect all negative function values
            forceToBePositive = []
            newGs = newGrid.getStorage()
            for i in xrange(newGs.size()):
                gp = newGs.get(i)
                if newNodalValues[newGs.seq(gp)] < 0.:
                    forceToBePositive.append(gp)

            if len(forceToBePositive) > 0:
                # if not, interpolate the function values for the negative
                # grid points
                if self.verbose:
                    warnings.warn("manually forcing the the function to be positive")
                newAlpha = self.makeCurrentGridPositive(newGrid, newAlpha)
                # newAlpha = SetGridPointsToZero().computeHierarchicalCoefficients(grid, newAlpha, forceToBePositive)
                # newAlpha = InterpolateParents().computeHierarchicalCoefficients(grid, newAlpha, forceToBePositive)

            if newGrid.getSize() == grid.getSize():
                break

            if self.verbose:
                fig = plt.figure()
                plotSG2d(newGrid, newAlpha)
                plt.title("iteration = %i" % iteration)
                fig.show()

                print "%i + %i = %i for discretization" % (grid.getSize(),
                                                           newGrid.getSize() - grid.getSize(),
                                                           newGrid.getSize(),)

            iteration += 1
            grid, alpha = newGrid, newAlpha

        # coarsening: remove all new grid points with negative surplus
        coarsedGrid, coarsedAlpha, _ = self.coarsening(newGrid, newAlpha, newGridPoints)
        if self.verbose:
            print "%i - %i = %i grid size after coarsening" % (newGrid.getSize(),
                                                               newGrid.getSize() - coarsedGrid.getSize(),
                                                               coarsedGrid.getSize(),)

        # scale the result such that the integral over it is one
#         c = createOperationQuadrature(coarsedGrid).doQuadrature(coarsedAlpha)
#         coarsedAlpha.mult(1. / c)

        # security check for positiveness
        checkPositivity(coarsedGrid, coarsedAlpha)

        return coarsedGrid, coarsedAlpha
