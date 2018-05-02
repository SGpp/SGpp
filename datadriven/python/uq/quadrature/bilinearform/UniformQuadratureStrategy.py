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

    def getKey(self, gps):
        """
        Generates a unique key for a given list of grid points
        @param gps: list of HashGridPoint
        """
        return tuple([(gp.getLevel(d), gp.getIndex(d)) for gp in gps for d in xrange(gp.getDimension())])

    def computeBilinearForm(self, grid):
        """
        Compute bilinear form for the current grid
        @param grid: Grid
        @return: DataMatrix
        """
        gs = grid.getStorage()
        A = DataMatrix(gs.getSize(), gs.getSize())
        A.setAll(0.)
        createOperationLTwoDotExplicit(A, grid)
        A = A.array()

        # store the result in the hash map
        for i in xrange(gs.getSize()):
            gpi = gs.getPoint(i)
            for j in xrange(gs.getSize()):
                gpj = gs.getPoint(j)
                key = self.getKey([gpi, gpj])
                self._map[key] = A[i, j]
        
        return A

    def computeBilinearFormEntry(self, basis, gpi, gpj):
        return self._map[self.getKey([gpi, gpj])], 0.
