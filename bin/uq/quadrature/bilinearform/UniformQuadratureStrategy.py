"""
Created on Aug 6, 2014

@author: franzefn
"""
from pysgpp import DataMatrix, createOperationLTwoDotExplicit
from BilinearQuadratureStrategy import BilinearQuadratureStrategy


class UniformQuadratureStrategy(BilinearQuadratureStrategy):
    """
    Generic object for quadrature strategies
    """

    def __init__(self):
        """
        Constructor
        """
        super(self.__class__, self).__init__()

    def computeBilinearForm(self, grid):
        """
        Compute bilinear form for the current grid
        @param grid: Grid
        @return: DataMatrix
        """
        gs = grid.getStorage()
        A = DataMatrix(gs.size(), gs.size())
        A.setAll(0.)
        createOperationLTwoDotExplicit(A, grid)
        for i in xrange(gs.size()):
            gpi = gs.get(i)
            for j in xrange(gs.size()):
                gpj = gs.get(j)
                key = self.getKey(gpi, gpj)
                self._map[key] = A.get(i, j)
        return A

    def computeBilinearFormEntry(self, basis, gpi, gpj):
        return self._map[self.getKey(gpi, gpj)], 0.
