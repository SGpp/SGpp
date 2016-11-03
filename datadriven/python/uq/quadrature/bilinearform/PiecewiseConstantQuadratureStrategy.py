"""
Created on Aug 6, 2014

@author: franzefn
"""
from pysgpp import DataMatrix, createOperationLTwoDotExplicit, DataVector
from BilinearQuadratureStrategy import BilinearQuadratureStrategy


class PiecewiseConstantQuadratureStrategy(BilinearQuadratureStrategy):
    """
    Generate the a quadrature strategy that appriximates the probability
    density in each dimension with a piecewise constant function. Each constant
    function is determined by the density evaluated at the center of the support
    of the corresponding basis function.
    """

    def __init__(self, params):
        """
        Constructor
        """
        super(self.__class__, self).__init__()
        self._U = params.getIndependentJointDistribution()

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
        @return DataMatrix
        """
        # create bilinear form of the grid
        gs = grid.getStorage()
        A = DataMatrix(gs.getSize(), gs.getSize())
        A.setAll(0.)
        createOperationLTwoDotExplicit(A, grid)

        # multiply the entries with the pdf at the center of the support
        p = DataVector(gs.getDimension())
        q = DataVector(gs.getDimension())

        for i in xrange(gs.getSize()):
            gpi = gs.getPoint(i)
            gs.getCoordinates(gpi, p)
            for j in xrange(gs.getSize()):
                gpj = gs.getPoint(j)
                gs.getCoordinates(gpj, q)
                y = float(A.get(i, j) * self._U.pdf(p))
                A.set(i, j, y)
                A.set(j, i, y)
                self._map[self.getKey([gpi, gpj])] = A.get(i, j)

        return A

    def computeBilinearFormEntry(self, basis, gpi, gpj):
        return self._map[self.getKey([gpi, gpj])], 0.
