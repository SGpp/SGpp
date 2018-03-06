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
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import parents,\
    getLevelIndex, getGridPointsOnBoundary


class ScaledMinOfParents(InterpolationAlgorithm):

    def computeMin(self, i, grid, alpha, nodalValues):
        gs = grid.getStorage()
        numDims = gs.getDimension()
        
        opEval = createOperationEval(grid)
        alphaVec = DataVector(alpha)
        
        gp = gs.getPoint(i)
        p = DataVector(numDims)
        gp.getStandardCoordinates(p)
        value = float("inf")
        
        level, index = getLevelIndex(gp)
        for idim in xrange(numDims):
            left, right = getGridPointsOnBoundary(level[idim], index[idim])
            
            if left is not None:
                llevel, lindex = left
                p[idim] = 2 ** -llevel * lindex
                leftValue = opEval.eval(alphaVec, p)
            else:
                leftValue = 0.0

            if right is not None:
                rlevel, rindex = right
                p[idim] = 2 ** -rlevel * rindex
                rightValue = opEval.eval(alphaVec, p)
            else:
                rightValue = 0.0

            interpolatedValue = abs(leftValue - rightValue) / 2

            if interpolatedValue < value:
                value = interpolatedValue

            # reset p
            p[idim] = 2 ** -level[idim] * index[idim]

        return value

    def computeHierarchicalCoefficients(self, grid, alpha, addedGridPoints):
        # define the order of computing the interpolated values -> this
        # makes sure that all function values of the hierarchical ancestors
        # exist for the current value we are interpolating
        nodalValues = dehierarchize(grid, alpha)
        neg = []
        for i, yi in enumerate(nodalValues):
            if yi < 0:
                neg.append(i)

        if len(neg) > 0:
            for i in neg:
                alpha[i] = alpha[i] - nodalValues[i] + self.computeMin(i, grid, alpha, nodalValues) - nodalValues[i]

#             alpha = hierarchize(grid, nodalValues)
#             print np.all(abs(newAlpha - alpha) < 1e-13)

            gs = grid.getStorage()
#             for i in xrange(len(alpha)):
#                 if i not in neg:
#                     if abs(alpha[i] - alphaNew[i]) > 1e-10:
#                         print "%i: %s -> |%g - %g| = %g" % (i,
#                                                             [gs.getPoint(i).getStandardCoordinate(d) for d in xrange(gs.getDimension())],
#                                                             alpha[i], alphaNew[i],
#                                                             abs(alpha[i] - alphaNew[i]))
# #                     assert abs(alpha[i] - alphaNew[i]) < 1e-10
#
#             alpha = alphaNew

            # check if the coefficients of the new grid points are positive
            if addedGridPoints is not None:
                for gp in addedGridPoints:
                    ix = gs.getSequenceNumber(gp)
                    if ix not in neg and alpha[ix] > 1e-13:
                        print "do not touch the non negative new grid points!!!"
                    if ix in neg and alpha[ix] < -1e-13:
                        print "negative coefficient found: %s -> %g (nodal=%g)" % ([gp.getStandardCoordinate(d) for d in xrange(gs.getDimension())],
                                                                                   alpha[ix],
                                                                                   nodalValues[ix])
#                     if ix in neg:
#                         assert alpha[ix] > -1e-13
        return alpha
