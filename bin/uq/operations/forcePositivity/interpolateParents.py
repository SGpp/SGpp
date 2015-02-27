'''
Created on Feb 6, 2015

@author: franzefn
'''

from interpolationAlgorithm import InterpolationAlgorithm
from pysgpp import createOperationEval, DataVector
from bin.uq.operations import (dehierarchize,
                               getBoundsOfSupport,
                               hierarchize)


class InterpolateParents(InterpolationAlgorithm):

    def computeMean(self, grid, alpha, newGridPoints):
        gs = grid.getStorage()
        opEval = createOperationEval(grid)

        newGs = grid.getStorage()
        newNodalValues = dehierarchize(grid, alpha)
        p = DataVector(gs.dim())

        for newGp in newGridPoints:
            # run over all grid points of current level
            newGp.getCoords(p)
            pp = DataVector(p)
            i = newGs.seq(newGp)
            newNodalValues[i] = 0.
            for d in xrange(gs.dim()):
                # get current index
                index = newGp.getIndex(d)
                level = newGp.getLevel(d)
                # bounds
                xlow, xhigh = getBoundsOfSupport(level, index)
                # compute function values at bounds
                pp[d] = xlow
                fxlow = opEval.eval(alpha, pp)
                pp[d] = xhigh
                fxhigh = opEval.eval(alpha, pp)
                # interpolate linearly
                a = (fxhigh - fxlow) / (xlow - xhigh)
                newNodalValues[i] += a * (p[d] - xlow) + fxhigh
                # reset pp and set sumWeights
                pp[d] = p[d]
            newNodalValues[i] /= gs.dim()

        return newNodalValues

    def computeHierarchicalCoefficients(self, grid, alpha, newGridPoints):
        # define the order of computing the interpolated values -> this
        # makes sure that all function values of the hierarchical ancestors
        # exist for the current value we are interpolating
        gs = grid.getStorage()
        levelsum = gs.getMaxLevel() + gs.dim()
        maxlevelsum = gs.getMaxLevel() * gs.dim()

        # interpolate all missing values
        while levelsum <= maxlevelsum:
            filteredGridPoints = []
            for newGp in newGridPoints:
                if levelsum == newGp.getLevelSum():
                    filteredGridPoints.append(newGp)

            # compute coefficients for new grid points
            if len(filteredGridPoints) > 0:
                nodalValues = self.computeMean(grid, alpha, filteredGridPoints)

            levelsum += 1

        return hierarchize(grid, nodalValues)
