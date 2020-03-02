# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

class InterpolationAlgorithm(object):

    def computeHierarchicalCoefficients(self, grid, alpha, newGridPoints):
        raise NotImplementedError()
