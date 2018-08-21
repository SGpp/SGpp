'''
Created on Feb 6, 2015

@author: franzefn
'''
from builtins import object


class InterpolationAlgorithm(object):

    def computeHierarchicalCoefficients(self, grid, alpha, newGridPoints):
        raise NotImplementedError()
