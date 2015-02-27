'''
Created on Feb 6, 2015

@author: franzefn
'''

from interpolationAlgorithm import InterpolationAlgorithm
from bin.uq.dists.SGDEdist import SGDEdist
import ConfigParser as cp
import numpy as np
from bin.uq.operations import dehierarchize
from uq.operations.forcePositivity.operationMakePositive import OperationMakePositive
from uq.operations.forcePositivity.interpolateParents import InterpolateParents


class EstimateDensityAlgorithm(InterpolationAlgorithm):

    def __init__(self, config):
        self.config = config

    def computeHierarchicalCoefficients(self, grid, alpha, newGridPoints):
        # load path and file names of the grid and if the coefficients from
        # the config file
        s = cp.ConfigParser()
        s.read(self.config)
        gridfile = s.get('dmest', 'writegridfile')
        alphafile = s.get('dmest', 'writealphafile')

        # serialize grid
        fd = open(gridfile, "w")
        fd.write(grid.serialize())
        fd.close()

        # serialize coefficients
        np.savetxt(alphafile, alpha.array())

        # estimate the density on the given grid
        dist = SGDEdist.byConfig(self.config)

        # return new new coefficient vector
        return dist.alpha
