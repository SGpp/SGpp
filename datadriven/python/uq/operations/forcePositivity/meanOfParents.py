'''
Created on Feb 6, 2015

@author: franzefn
'''

import numpy as np

from interpolationAlgorithm import InterpolationAlgorithm
from pysgpp import createOperationEval, DataVector
from pysgpp.extensions.datadriven.uq.operations import (dehierarchize,
                               getBoundsOfSupport,
                               hierarchize)
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import parents


class MeanOfParents(InterpolationAlgorithm):

    def computeMean(self, i, grid, nodalValues):
        gs = grid.getStorage()

        value = 0.0
        gpps = parents(grid, gs.getPoint(i))
        for _, gpp in gpps:
            value += nodalValues[gs.getSequenceNumber(gpp)]
        
        return 1. / len(gpps) * value

    def computeHierarchicalCoefficients(self, grid, alpha, addedGridPoints):
        # define the order of computing the interpolated values -> this
        # makes sure that all function values of the hierarchical ancestors
        # exist for the current value we are interpolating
        nodalValues = dehierarchize(grid, alpha)
        neg = np.array([], dtype="int")
        for i, yi in enumerate(nodalValues):
            if yi < 0:
                nodalValues[i] = 0
                neg = np.append(neg, i)

        if len(neg) > 0:
            for i in neg:
                # apply factor of 1/4 -> grid is in the convergent phase
                nodalValues[i] = 0.25 * self.computeMean(i, grid, nodalValues)

            alpha = hierarchize(grid, nodalValues)

            # check if the coefficients of the new grid points are positive
            if addedGridPoints is not None:
                gs = grid.getStorage()
                assert all([alpha[gs.getSequenceNumber(gp)] > -1e-13 for gp in addedGridPoints])

        return alpha
