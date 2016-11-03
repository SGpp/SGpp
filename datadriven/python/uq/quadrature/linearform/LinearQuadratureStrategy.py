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
        key = self._map.getKey([gpi], d)
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
            # run over all dimensions
            for d in xrange(gpi.getDimension()):
                # compute linear form for one entry
                value, erri = self.getLinearFormEntry(gs, gpi, basis, d)
                # collect results
                b[i] *= value
                err += erri
        return b, err

    def getLinearFormEntry(self, gs, gp, basis, d):
        """
        Restore the bilinear form of two grid points if it is available.
        If not, forward the result to the computation method.
        @param gs: HashGridStorage
        @param gp: HashGridPoint
        @param basis: SG++ Basis
        @param d: int dimension
        """
        available, key = self.hasValue(gp, d)
        if not available:
            # there is no information available for the current combination
            # of grid points
            val, err = self.computeLinearFormEntry(gs, gp, basis, d)
            # store value
            self._map[key] = val, err
        else:
            val, err = self._map[key]
        return val, err

    def computeLinearFormEntry(self, gs, gp, basis, d):
        """
        Compute the bilinear form of one grid point with another one
        @param gs: HashGridStorage
        @param gp: HashGridPoint
        @param basis: SG++ Basis
        """
        raise NotImplementedError()
