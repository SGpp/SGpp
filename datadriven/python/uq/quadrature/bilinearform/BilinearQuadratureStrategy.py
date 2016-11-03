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
        key = self._map.getKey([gpi, gpj], d)
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
        @param gps: list of HashGridPoint
        @param basisj: SG++ Basis for grid indices j
        @return DataVector
        """
        b = np.ones(len(gpsj))
        err = 0.
        # run over all entries
        for j, gpj in enumerate(gpsj):
            # run over all dimensions
            for d in xrange(gpi.getDimension()):
                # compute bilinear form for one entry
                s, erri = self.getBilinearFormEntry(gs, gpi, basisi, gpj, basisj, d)
                # combine different dimensions
                b[j] *= s
                err += erri
        return b, err

    def getBilinearFormEntry(self, gs, gpi, basisi, gpj, basisj, d):
        """
        Restore the bilinear form of two grid points if it is available.
        If not, forward the result to the computation method.
        @param gs: HashGridStorage
        @param gpi: HashGridPoint
        @param basisi: SG++ Basis
        @param gpj: HashGridPoint
        @param basisj: SG++ Basis
        @param d: int dimension
        """
        available, key = self.hasValue(gpi, gpj, d)
        if not available:
            # there is no information available for the current combination
            # of grid points
            val, err = self.computeBilinearFormEntry(gs, gpi, basisi, gpj, basisj, d)
            # store value
            self._map[key] = val, err
        else:
            val, err = self._map[key]
        return val, err

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
