# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.operations.forcePositivity.interpolationAlgorithm import InterpolationAlgorithm
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
