#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org
#
from pysgpp import Grid, DataVector
from bin.uq.transformation import LinearTransformation
from bin.uq.operations import hierarchize, discretize
from bin.uq.operations.discretization import discretizeFunction
"""
@file    base.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 14:04:16 2013

@brief   Superclass for all distributions

@version  0.1
"""


class Dist(object):
    """
    The Dist class, which is the super class for all
    distributions of this package
    """

    def pdf(self, p, *args, **kws):
        """
        Probability distribution function
        @param p: (tuple) of floats
        @return: probability distribution value
        """
        raise NotImplementedError()

    def cdf(self, p, *args, **kws):
        """
        Cumulative distribution function
        @param p: (tuple) of floats
        @return: cumulative distribution value
        """
        raise NotImplementedError()

    def ppf(self, p, *args, **kws):
        """
        Point percentile function
        @param p: (tuple) of floats
        @return: point percentile value
        """
        raise NotImplementedError()

    def rvs(self, n=1):
        """
        Generates n random numbers w.r.t. the marginal distributions
        @param n: int number of random values
        @return: numpy array [n, dim]
        """
        raise NotImplementedError()

    def mean(self):
        """
        @return: expectation value
        """
        raise NotImplementedError()

    def var(self):
        """
        @return: variance
        """
        raise NotImplementedError()

    def std(self):
        """
        @return: standard deviation
        """
        raise NotImplementedError()

    def getBounds(self):
        """
        Get the distribution's intervals
        @return: numpy array [dim, 2]
        """
        raise NotImplementedError()

    def getDim(self):
        """
        Get number of marginal distributions
        @return: int number of marginal distributions
        """
        raise NotImplementedError()

    def discretize(self, *args, **kws):
        """
        discretize the pdf of the current distribution
        using a sparse grid interpolant
        """
        bounds = self.getBounds()
        if self.getDim() == 1:
            bounds = [bounds]
        return discretizeFunction(self.pdf, bounds, hasBorder=False,
                                  *args, **kws)

    @classmethod
    def fromJson(cls, jsonObject):
        import bin.uq.dists as dists
        if jsonObject['module'] == 'bin.uq.dists.J':
            return dists.J.fromJson(jsonObject)
        elif jsonObject['module'] == 'bin.uq.dists.Corr':
            return dists.Corr.fromJson(jsonObject)
        elif jsonObject['module'] == 'bin.uq.dists.SGDEdist':
            return dists.SGDEdist.fromJson(jsonObject)
        elif jsonObject['module'] == 'bin.uq.dists.Uniform':
            return dists.Uniform.fromJson(jsonObject)
        elif jsonObject['module'] == 'bin.uq.dists.TNormal':
            return dists.TNormal.fromJson(jsonObject)
        elif jsonObject['module'] == 'bin.uq.dists.Normal':
            return dists.Normal.fromJson(jsonObject)
        elif jsonObject['module'] == 'bin.uq.dists.Lognormal':
            return dists.Lognormal.fromJson(jsonObject)
        elif jsonObject['module'] == 'bin.uq.dists.Beta':
            return dists.Beta.fromJson(jsonObject)
        elif jsonObject['module'] == 'bin.uq.dists.MultivariateNormal':
            return dists.MultivariateNormal.fromJson(jsonObject)
        else:
            raise TypeError('Unknown distribution "%s" => Please register \
                             it in fromJson function' % jsonObject['module'])
