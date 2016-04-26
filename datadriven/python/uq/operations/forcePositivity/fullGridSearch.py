from pysgpp import Grid, DataVector, createOperationEval, HashGridIndex
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

        opEval = createOperationEval(grid)
        alphaVec = DataVector(alpha)
        p = DataVector(self.numDims)
        for i in xrange(fullGridStorage.getSize()):
            gp = fullGridStorage.get(i)
            gp.getCoords(p)
            if not gs.has_key(gp) and opEval.eval(alphaVec, p) < 0.0:
                self.candidates.append(HashGridIndex(gp))
