"""
Created on Aug 6, 2014

@author: franzefn
"""
from pysgpp import DataMatrix, DataVector
from pysgpp.extensions.datadriven.uq.quadrature.HashQuadrature import HashQuadrature


class TrilinearQuadratureStrategy(HashQuadrature):
    """
    Generic object for quadrature strategies
    """
    def hasValue(self, gpk, gpi, gpj, d):
        key = self._map.getKey([gpk, gpi, gpj], d)
        if key in self._map:
            return True, key
        key = self._map.getKey([gpk, gpj, gpi], d)
        if key in self._map:
            return True, key
        return False, key

    def computeTrilinearFormByList(self,
                                   gpsk, basisk, alphak,
                                   gpsi, basisi,
                                   gpsj, basisj):
        """
        Compute trilinear form for two lists of grid points
        @param gpsk: list of HashGridIndex
        @param basisk: SG++ basis for grid indices gpsk
        @param alphak: coefficients for kth grid
        @param gpsi: list of HashGridIndex
        @param basisi: SG++ basis for grid indices gpsi
        @param gpsj: list of HashGridIndex
        @param basisj: SG++ basis for grid indices gpsj
        @return: DataMatrix
        """
        print "# evals: %i^2 * %i = %i" % (len(gpsi), len(gpsk), len(gpsi) ** 2 * len(gpsk))
        A = DataMatrix(len(gpsi), len(gpsj))
        err = 0.
        # run over all rows
        for i, gpi in enumerate(gpsi):
            # run over all columns
            for j, gpj in enumerate(gpsj):
                # run over all gpks
                b, erri = self.computeTrilinearFormByRow(gpsk, basisk,
                                                         gpi, basisi,
                                                         gpj, basisj)
                # get the overall contribution in the current dimension
                value = alphak.dotProduct(b)
                A.set(i, j, value)

                # error statistics
                err += erri
        return A, err

    def computeTrilinearFormByRow(self,
                                  gpsk, basisk,
                                  gpi, basisi,
                                  gpj, basisj):
        """
        Compute the trilinear form of two grid point with a list
        of grid points
        @param gpk: list of HashGridIndex
        @param basisk: SG++ Basis for grid indices k
        @param gpi: HashGridIndex
        @param basisi: SG++ Basis for grid indices i
        @param gpj: HashGridIndex
        @param basisj: SG++ Basis for grid indices j
        @return DataVector
        """
        b = DataVector(len(gpsk))
        b.setAll(1.0)
        err = 0.
        # run over all entries
        for k, gpk in enumerate(gpsk):
            # run over all dimensions
            for d in xrange(gpi.getDimension()):
                # compute trilinear form for one entry
                value, erri = self.getTrilinearFormEntry(gpk, basisk,
                                                         gpi, basisi,
                                                         gpj, basisj,
                                                         d)
                b[k] *= value
                err += erri
        return b, err

    def getTrilinearFormEntry(self,
                              gpk, basisk,
                              gpi, basisi,
                              gpj, basisj,
                              d):
        """
        Restore the trilinear form of two grid points if it is available.
        If not, forward the result to the computation method.
        @param gpk: HashGridIndex
        @param basisk: SG++ Basis
        @param gpi: HashGridIndex
        @param basisi: SG++ Basis
        @param gpj: HashGridIndex
        @param basisj: SG++ Basis
        """
        available, key = self.hasValue(gpk, gpi, gpj, d)
        if not available:
            # there is no information available for the current combination
            # of grid points
            val, err = self.computeTrilinearFormEntry(gpk, basisk,
                                                      gpi, basisi,
                                                      gpj, basisj,
                                                      d)
            # store value
            self._map[key] = val, err
        else:
            val, err = self._map[key]
        return val, err

    def computeTrilinearFormEntry(self, gpk, basisk, gpi, basisi, gpj, basisj, d):
        """
        Compute the Trilinear form of one grid point with another one
        @param gpk: HashGridIndex
        @param basisk: SG++ Basis
        @param gpi: HashGridIndex
        @param basisi: SG++ Basis
        @param gpj: HashGridIndex
        @param basisj: SG++ Basis
        """
        raise NotImplementedError()
