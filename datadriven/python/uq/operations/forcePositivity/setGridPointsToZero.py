'''
Created on Feb 6, 2015

@author: franzefn
'''

from interpolationAlgorithm import InterpolationAlgorithm
from pysgpp import createOperationEval, DataVector
from pysgpp.extensions.datadriven.uq.operations import (dehierarchize,
                               getBoundsOfSupport,
                               hierarchize)


class SetGridPointsToZero(InterpolationAlgorithm):

    def computeHierarchicalCoefficients(self, grid, alpha, newGridPoints):
        newGs = grid.getStorage()
        nodalValues = dehierarchize(grid, alpha)

        for newGp in newGridPoints:
            # run over all grid points of current level
            i = newGs.getSequenceNumber(newGp)
            nodalValues[i] = 0.

        return hierarchize(grid, nodalValues)
