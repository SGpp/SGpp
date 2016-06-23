
'''
Created on Feb 6, 2015

@author: franzefn
'''

import os
import numpy as np

from interpolationAlgorithm import InterpolationAlgorithm
from pysgpp.extensions.datadriven.uq.dists.SGDEdist import SGDEdist
import tempfile, uuid, json
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBasis, \
    dehierarchize, getLevelIndex, hierarchize
from pysgpp import DataVector,createOperationLTwoDotExplicit, RegularizationType_Identity, RegularizationType_Laplace
from pysgpp._pysgpp_swig import createOperationIdentity, createOperationLaplace


class EstimateDensityAlgorithm(InterpolationAlgorithm):

    def __init__(self, trainSamples, lmbd, regularizationType):
        self.trainSamples = trainSamples
        self.lmbd = lmbd
        self.regularizationType = regularizationType

    def computeRegularizationMatrix(self, grid):
        if self.regularizationType == RegularizationType_Identity:
            return createOperationIdentity(grid)
        elif self.regularizationType == RegularizationType_Laplace:
            return createOperationLaplace(grid)
        else:
            raise AttributeError("EstiamateDensityAlgorithm: regularization type is unknown")

    
    def computeL2Matrix(self, grid):
        return createOperationLTwoDotExplicit(grid)


    def evalBasis(self, level, index, basis, x):
        return np.prod([basis.eval(level[idim], index[idim], x[idim])
                        for idim in xrange(len(x))])
    
    def evalL2(self, A, alpha):
        result = DataVector(len(alpha))
        alpha_vec = DataVector(alpha)
        A.mult(alpha_vec, result)
        return result.array()

    def evalRegularization(self, C, alpha):
        result = DataVector(len(alpha))
        alpha_vec = DataVector(alpha)
        C.mult(alpha_vec, result)
        return result.array()
    
    def evaluateDensityLocally(self, ix, grid, A, C, a, c):
        # prepare variables
        gs = grid.getStorage()
        numDims = gs.getDimension()

        gp = gs.get(ix)
        levels, indices = getLevelIndex(gp)
        basis = grid.getBasis()

        # compute right hand side
        rhs = 0
        for i in xrange(self.trainSamples.shape[0]):
            rhs += self.evalBasis(levels, indices, basis, self.trainSamples[i, :])
        rhs /= 1.0 * self.trainSamples.shape[0]

        # compute left hand side A alpha + lambda C alpha
        lhs = a[ix] + self.lmbd * c[ix]

        if rhs - lhs < 0.0:
            return 0.0
        else:
            alpha = np.zeros(gs.getSize())
            alpha[ix] = 1.0
            aix = self.evalL2(A, alpha)[ix]
            cix = self.evalRegularization(C, alpha)[ix]

            return (rhs - lhs)  # / (aix + self.lmbd * cix)


    def computeHierarchicalCoefficients(self, grid, alpha, addedGridPoints=None):
        nodalValues = dehierarchize(grid, alpha)
        neg = []
        for i, yi in enumerate(nodalValues):
            if yi < 0:
                neg.append(i)
                nodalValues[i] = 0.0

        if len(neg) > 0:
            # compute the coefficients for each grid point by estimating
            # the local density
#             C = self.computeRegularizationMatrix(grid)
#             A = self.computeL2Matrix(grid)
#             a = self.evalL2(A, alpha)
#             c = self.evalRegularization(C, alpha)

            for i in neg:
#                 nodalValues[i] = self.evaluateDensityLocally(i, grid, A, C, a, c)
                # compute right hand side
                for i in xrange(self.trainSamples.shape[0]):
                    nodalValues[i] += self.evalBasis(levels, indices, basis, self.trainSamples[i, :])
                nodalValues[i] /= 1.0 * self.trainSamples.shape[0]

            alpha = hierarchize(grid, nodalValues)

            # check if the coefficients of the new grid points are positive
            if addedGridPoints is not None:
                gs = grid.getStorage()
                assert all([alpha[gs.getSequenceNumber(gp)] > -1e-13 for gp in addedGridPoints])
        return alpha
