"""
Created on Aug 6, 2014

@author: franzefn
"""
from pysgpp.extensions.datadriven.uq.quadrature.HashQuadrature import HashQuadrature
import numpy as np


class BilinearQuadratureStrategy(HashQuadrature):
    """
    Generic object for quadrature strategies
    """

    def hasValue(self, gpi, gpj, d):
        key = self._map.getKey([gpi, gpj], d)
        if key in self._map:
            return True, key
        return False, key

    def computeBilinearFormByList(self, gpsi, basisi, gpsj, basisj):
        """
        Compute bilinear form for two lists of grid points
        @param gpsi: list of HashGridIndex
        @param basisi: SG++ basis for grid indices gpsi
        @param gpsj: list of HashGridIndex
        @param basisj: SG++ basis for grid indices gpsj
        @return: numpy array
        """
        A = np.ndarray((len(gpsi), len(gpsj)))
        err = 0.
        # run over all rows
        for i, gpi in enumerate(gpsi):
            b, erri = self.computeBilinearFormByRow(gpi, basisi, gpsj, basisj)
            A[i, :] = b
            err += erri
        return A, err

    def computeBilinearFormByRow(self, gpi, basisi, gpsj, basisj):
        """
        Compute the bilinear form of one grid point with a list
        of grid points
        @param gpi: HashGridIndex
        @param basisi: SG++ Basis for grid indices i
        @param gps: list of HashGridIndex
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
                s, erri = self.getBilinearFormEntry(gpi, basisi, gpj, basisj, d)
                # combine different dimensions
                b[j] *= s
                err += erri
        return b, err

    def getBilinearFormEntry(self, gpi, basisi, gpj, basisj, d):
        """
        Restore the bilinear form of two grid points if it is available.
        If not, forward the result to the computation method.
        @param gpi: HashGridIndex
        @param basisi: SG++ Basis
        @param gpj: HashGridIndex
        @param basisj: SG++ Basis
        @param d: int dimension
        """
        available, key = self.hasValue(gpi, gpj, d)
        if not available:
            # there is no information available for the current combination
            # of grid points
            val, err = self.computeBilinearFormEntry(gpi, basisi, gpj, basisj, d)
            # store value
            self._map[key] = val, err
        else:
            val, err = self._map[key]
        return val, err

    def computeBilinearFormEntry(self, gpi, basisi, gpj, basisj, d):
        """
        Compute the bilinear form of one grid point with another one
        @param gpi: HashGridIndex
        @param basisi: SG++ Basis
        @param gpj: HashGridIndex
        @param basisj: SG++ Basis
        @param d: int dimension
        """
        raise NotImplementedError()
