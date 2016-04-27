from pysgpp import Grid, DataVector, createOperationEval, HashGridIndex
from findCandidateSet import CandidateSet
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import parent, \
    hasAllChildren, getHierarchicalAncestors


class ChildrenNodesAsCandidates(CandidateSet):


    def __init__(self, grid, alpha):
        super(ChildrenNodesAsCandidates, self).__init__()


    def findCandidates(self, grid, alpha, addedGridPoints):
        """
        collect all leaf nodes with at least one negative hierarchical
        ancestor
        """
        # collect all leaf nodes with at least one negative hierarchical
        # ancestor
        refinementCandidates = {}
        gs = grid.getStorage()
        numDims, numGridPoints = gs.getDimension(), gs.getSize()

        self.costs = 0
        if self.iteration == 0:
            for i in xrange(numGridPoints):
                self.costs += 1
                gp = gs.get(i)
                if not hasAllChildren(grid, gp):
                    refinementCandidates[i] = gp
        else:
            for gp in self.candidates:
                self.costs += 1
                if not hasAllChildren(grid, gp):
                    refinementCandidates[i] = gp

        # add all children to the list of candidates
        maxLevel = gs.getMaxLevel()
        for i, gp in refinementCandidates.items():
            self.costs += 1
            for idim in xrange(numDims):
                if gp.getLevel(idim) < maxLevel:
                    gpl = HashGridIndex(gp)
                    gs.left_child(gpl, idim)
                    if not gs.has_key(gpl):
                        self.candidates.append(gpl)

                    # get right child
                    gpr = HashGridIndex(gp)
                    gs.right_child(gpr, idim)
                    if not gs.has_key(gpr):
                        self.candidates.append(gpr)
