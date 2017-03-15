"""
Created on Aug 6, 2014

@author: franzefn
"""
from pysgpp import DataMatrix, createOperationLTwoDotExplicit
from BilinearQuadratureStrategy import BilinearQuadratureStrategy


class UniformQuadratureStrategy(BilinearQuadratureStrategy):
    """
    Generate the a quadrature strategy for uniformly distributed
    random variables. Therefore is the mass matrix just the bilinear form
    of the sparse grid basis functions.
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

        # store the result in the hash map
        for i in xrange(gs.size()):
            gpi = gs.getPoint(i)
            for j in xrange(gs.size()):
                gpj = gs.getPoint(j)
                key = self.getKey(gpi, gpj)
                self._map[key] = A.get(i, j)
        return A

    def computeBilinearFormEntry(self, basis, gpi, gpj):
        return self._map[self.getKey(gpi, gpj)], 0.
