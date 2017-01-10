"""
Created on Aug 6, 2014

@author: franzefn
"""
from pysgpp.extensions.datadriven.uq.operations import getBasis
from pysgpp.extensions.datadriven.uq.quadrature import HashQuadrature
import numpy as np


class LinearQuadratureStrategy(HashQuadrature):
    """
    Generic object for quadrature strategies
    """

    def hasValue(self, gpi, d):
        if len(self._U) < d:
            U = self._U[d]
        else:
            U = self._U[-1]

        key = self._map.getKey(U, [gpi], d)
        if key in self._map:
            return True, key
        return False, key

    def computeLinearForm(self, grid):
        """
        Compute bilinear form for the current grid
        @param grid: Grid
        @return numpy array
        """
        gs = grid.getStorage()
        basis = getBasis(grid)
        v = np.ndarray(gs.size())
        err = 0.
        # run over all rows
        for i in xrange(gs.size()):
            gpi = gs.getPoint(i)
            # compute bilinear form for one entry
            v[i], erri = self.getLinearFormEntry(gs, gpi, basis)
            err += erri
        return v, err

    def computeLinearFormByList(self, gs, gps, basis):
        """
        Compute bilinear form for two lists of grid points
        @param gs: HashGridStorage
        @param gps: list of HashGridPoint
        @param basis: SG++ basis for grid indices gpsi
        @return: numpy array
        """
        b = np.ones(len(gps))
        err = 0.
        # run over all items
        for i, gpi in enumerate(gps):
            # compute linear form for one entry
            value, erri = self.getLinearFormEntry(gs, gpi, basis)

            # collect results
            b[i] = value
            err += erri
        return b, err

    def getLinearFormEntry(self, gs, gp, basis):
        """
        Restore the bilinear form of two grid points if it is available.
        If not, forward the result to the computation method.
        @param gs: HashGridStorage
        @param gp: HashGridPoint
        @param basis: SG++ Basis
        """
        ans, err = 1.0, 0.0

        # run over all dimensions
        for d in xrange(gp.getDimension()):
            # compute linear form for one entry
            available, keyd = self.hasValue(gp, d)
            if not available:
                val, erri = self.computeLinearFormEntry(gs, gp, basis, d)
                # store value
                self._map[keyd] = val, erri
            else:
                val, erri = self._map[keyd]

            # collect results
            ans *= val
            err += erri

        return ans, err

    def computeLinearFormEntry(self, gs, gp, basis, d):
        """
        Compute the bilinear form of one grid point with another one
        @param gs: HashGridStorage
        @param gp: HashGridPoint
        @param basis: SG++ Basis
        """
        raise NotImplementedError()
