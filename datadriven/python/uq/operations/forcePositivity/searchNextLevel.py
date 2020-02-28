# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp import Grid, DataVector, createOperationEval, HashGridPoint
from pysgpp.extensions.datadriven.uq.operations.forcePositivity.findCandidateSet import CandidateSet
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import parent, \
    hasAllChildren, getHierarchicalAncestors, getLevel, getIndex


class SearchLevelWiseForCandidates(CandidateSet):


    def __init__(self, grid, alpha):
        super(SearchLevelWiseForCandidates, self).__init__()
        self.minLevelSum = 1e10
        self.maxLevelSum = -1


    def getAllChildrenNodesUpToMaxLevel(self, gp, maxLevel, grid):
        children = {}
        gs = grid.getStorage()
        
        for idim in range(gp.getDimension()):
            if gp.getLevel(idim) < maxLevel:
                self.costs += 1
                gpl = HashGridPoint(gp)
                gp.getLeftChild(idim)
                level, index = tuple(getLevel(gpl)), tuple(getIndex(gpl))
                children[level, index] = gpl

                # get right child
                self.costs += 1
                gpr = HashGridPoint(gp)
                gpr.getRightChild(idim)
                level, index = tuple(getLevel(gpr)), tuple(getIndex(gpr))
                children[level, index] = gpr
        return children
    

    def findCandidates(self, grid, alpha, addedGridPoints):
        """
        collect all leaf nodes with at least one negative hierarchical
        ancestor
        """
        # collect all leaf nodes with at least one negative hierarchical
        # ancestor
        refinementCandidates = []
        gs = grid.getStorage()
        numDims, numGridPoints = gs.getDimension(), gs.getSize()
        maxLevel = gs.getMaxLevel()

        self.costs = 0
        if self.iteration == 0:
            candidates = []
            for i in range(numGridPoints):
                self.costs += 1
                gp = gs.getPoint(i)
                if alpha[i] < 0.0:
                    candidates.append(HashGridPoint(gp))
                    self.minLevelSum = min(self.minLevelSum, gp.getLevelSum())
                    self.maxLevelSum = max(self.maxLevelSum, gp.getLevelSum())

            # refine all the grid points up to the maximum level sum
            # -> let the search progressing uniformly such that
            #    we look at grid points with the same level sum at most once
            #    to achieve the minimum amount of grid points
            # this corresponds basically to the search of all leaf nodes with
            # an ancestor of negative coefficient
            refinementCandidates = []
            for gp in candidates:
                diff = self.maxLevelSum - gp.getLevelSum()
                if diff == 0:
                    refinementCandidates.append(gp)
                else:
                    newCandidates = [gp]
                    for i in range(diff):
                        while len(newCandidates) > 0:
                            candidate = newCandidates.pop()
                            children = self.getAllChildrenNodesUpToMaxLevel(candidate, maxLevel, grid)
                            for childCandidate in list(children.values()):
                                if childCandidate.getLevelSum() == self.maxLevelSum:
                                    refinementCandidates.append(childCandidate)
                                elif childCandidate.getLevelSum() < self.maxLevelSum:
                                    newCandidates.append(childCandidate)
        else:
            refinementCandidates = list(self.newCandidates.values())

        self.newCandidates = {}
        self.candidates = []
        for gp in refinementCandidates:
            children = self.getAllChildrenNodesUpToMaxLevel(gp, maxLevel, grid)
            for (level, index), ngp in list(children.items()):
                if (level, index) not in self.newCandidates:
                    if not gs.isContaining(ngp):
                        self.candidates.append(ngp)
                    self.newCandidates[level, index] = ngp
