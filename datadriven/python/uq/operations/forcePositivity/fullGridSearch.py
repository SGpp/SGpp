from pysgpp import Grid, DataVector, createOperationEval, HashGridPoint
from findCandidateSet import CandidateSet
import numpy as np


class FullGridCandidates(CandidateSet):


    def __init__(self, grid):
        super(FullGridCandidates, self).__init__()
        # genreate new full grid
        gs = grid.getStorage()
        maxLevel = gs.getMaxLevel()
        self.numDims = gs.getDimension()

        self.fullGrid = Grid.createLinearGrid(self.numDims)
        self.fullGrid.getGenerator().full(maxLevel)


    def findCandidates(self, grid, alpha, addedGridPoints):
        fullGridStorage = self.fullGrid.getStorage()
        gs = grid.getStorage()

        if self.iteration == 0:
            self.costs += fullGridStorage.getSize()
        elif len(addedGridPoints) == 0:
            return
        
        opEval = createOperationEval(grid)
        for i in xrange(fullGridStorage.getSize()):
            gp = fullGridStorage.getPoint(i)
            if not gs.isContaining(gp):
                self.candidates.append(HashGridPoint(gp))
