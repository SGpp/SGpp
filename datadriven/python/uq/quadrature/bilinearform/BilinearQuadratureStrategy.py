"""
Created on Aug 6, 2014

@author: franzefn
"""
from pysgpp.extensions.datadriven.uq.quadrature.HashQuadrature import HashQuadrature
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBasis


class BilinearQuadratureStrategy(HashQuadrature):
    """
    Generic object for quadrature strategies
    """

    def hasValue(self, gpi, gpj, d):
        if len(self._U) < d:
            U = self._U[d]
        else:
            U = self._U[-1]

        key = self._map.getKey(U, [gpi, gpj], d)
        if key in self._map:
            return True, key
        return False, key

    def computeBilinearForm(self, grid):
        """
        Compute bilinear form for the current grid
        @param grid: Grid
        @return DataMatrix
        """
        # create bilinear form of the grid
        gs = grid.getStorage()

        gps = [None] * gs.getSize()
        for i in xrange(gs.getSize()):
            gps[i] = gs.getPoint(i)

        basis = getBasis(grid)
        A, err = self.computeBilinearFormByList(gs, gps, basis, gps, basis)

        return A, err

    def computeBilinearFormByList(self, gs, gpsi, basisi, gpsj, basisj):
        """
        Compute bilinear form for two lists of grid points
        @param gs: HashGridStorage
        @param gpsi: list of HashGridPoint
        @param basisi: SG++ basis for grid indices gpsi
        @param gpsj: list of HashGridPoint
        @param basisj: SG++ basis for grid indices gpsj
        @return: numpy array
        """
        A = np.ndarray((len(gpsi), len(gpsj)))
        err = 0.
        # run over all rows
        for i, gpi in enumerate(gpsi):
            b, erri = self.computeBilinearFormByRow(gs, gpi, basisi, gpsj, basisj)
            A[i, :] = b
            err += erri
        return A, err

    def computeBilinearFormByRow(self, gs, gpi, basisi, gpsj, basisj):
        """
        Compute the bilinear form of one grid point with a list
        of grid points
        @param gs: HashGridStorage
        @param gpi: HashGridPoint
        @param basisi: SG++ Basis for grid indices i
        @param gpsj: list of HashGridPoint
        @param basisj: SG++ Basis for grid indices j
        @return DataVector
        """
        b = np.ones(len(gpsj))
        err = 0.
        # run over all items
        for j, gpj in enumerate(gpsj):
            # compute bilinear form for one entry
            value, erri = self.getBilinearFormEntry(gs, gpi, basisi, gpj, basisj)

            # collect results
            b[j] = value
            err += erri

        return b, err


    def getBilinearFormEntry(self, gs, gpi, basisi, gpj, basisj):
        """
        Restore the bilinear form of two grid points if it is available.
        If not, forward the result to the computation method.
        @param gs: HashGridStorage
        @param gpi: HashGridPoint
        @param basisi: SG++ Basis
        @param gpj: HashGridPoint
        @param basisj: SG++ Basis
        """
        ans, err = 1.0, 0.0

        # run over all dimensions
        for d in xrange(gpi.getDimension()):
            # compute linear form for one entry
            available, keyd = self.hasValue(gpi, gpj, d)
            if not available:
                val, erri = self.computeBilinearFormEntry(gs, gpi, basisi, gpj, basisj, d)
                # store value
                self._map[keyd] = val, erri
            else:
                val, erri = self._map[keyd]

            # collect results
            ans *= val
            err += erri

        return ans, err


    def computeBilinearFormEntry(self, gs, gpi, basisi, gpj, basisj, d):
        """
        Compute the bilinear form of one grid point with another one
        @param gs: grid storage
        @param gpi: HashGridPoint
        @param basisi: SG++ Basis
        @param gpj: HashGridPoint
        @param basisj: SG++ Basis
        @param d: int dimension
        """
        raise NotImplementedError()
