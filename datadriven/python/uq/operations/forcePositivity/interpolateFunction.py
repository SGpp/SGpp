# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.operations.forcePositivity.interpolationAlgorithm import InterpolationAlgorithm
from pysgpp import createOperationEval, DataVector
from pysgpp.extensions.datadriven.uq.operations import (dehierarchize,
                                                        hierarchize)
import numpy as np


class InterpolateFunction(InterpolationAlgorithm):

    def __init__(self, func):
        self.func = func

    def computeHierarchicalCoefficients(self, grid, alpha, newGridPoints):
        # define the order of computing the interpolated values -> this
        # makes sure that all function values of the hierarchical ancestors
        # exist for the current value we are interpolating
        gs = grid.getStorage()
        nodalValues = np.ndarray(gs.getSize())
        p = DataVector(gs.getDimension())
        for i in range(gs.getSize()):
            gs.getPoint(i).getStandardCoordinates(p)
            nodalValues[i] = self.func(p.array())

        return hierarchize(grid, nodalValues)
