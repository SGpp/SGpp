#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    base.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 14:04:16 2013

@brief   Superclass for all distributions

@version  0.1
"""

from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation
from pysgpp.extensions.datadriven.uq.sampler.Sample import Sample, SampleType
from pysgpp import DataVector

import numpy as np


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
        return np.sqrt(self.var())

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

    def cov(self):
        """
        Get covariance matrix
        """
        return np.diag(np.ones(self.getDim()) * self.var())

    def corrcoeff(self, covMatrix=None):
        """
        """
        if covMatrix is None:
            covMatrix = self.cov()

        numDims = covMatrix.shape[0]
        corr = np.ndarray(covMatrix.shape)
        for idim in xrange(numDims):
            sigmai = np.sqrt(covMatrix[idim, idim])
            for jdim in xrange(idim + 1, numDims):
                sigmaj = np.sqrt(covMatrix[jdim, jdim])
                corrij = covMatrix[idim, jdim] / (sigmai * sigmaj)
                corr[idim, jdim] = corr[jdim, idim] = corrij
            corr[idim, idim] = 1.0
        return corr

    def klDivergence(self, dist, testSamplesUnit=None, testSamplesProb=None,
                     n=1e4):
        r"""
        computes the KL-divergence from this distribution with respect to dist

        \approx \frac{1}{n} \sum_{i = 1}^{n} p(x_i) log_2 (p(x_i) / q(x_i))

        and for samples obtained via importance sampling it holds

        \approx \frac{1}{n} \sum_{i = 1}^{n} log_2 (p(y_i) / q(y_i))
        = \frac{1}{n} \sum_{i = 1}^{n} log_2 p(y_i) - log_2 q(y_i))
        = [\frac{1}{n} \sum_{i = 1}^{n} log_2 p(y_i)] - [\frac{1}{n} \sum_{i = 1}^{n} log_2 q(y_i)]
        = mean(log_2 p(y_i)) - mean(log_2 q(y_i))

        @param dist: Dist
        @param testSamplesUnit: numpy array
        @param testSamplesProb: numpy array
        """
        # compute log_2(p(x_i)) and log_2(q(x_i))
        n = testSamplesUnit.shape[0]

        if testSamplesProb is None:
            testSamplesProb = testSamplesUnit

        p = np.zeros(n)
        q = np.zeros(n)
        for i in xrange(n):
            p[i] = self.pdf(testSamplesProb[i, :])
            q[i] = max(1e-10, dist.pdf(testSamplesUnit[i, :]))

        # compute the KL divergence
        ans = np.mean(np.log2(p)) - np.mean(np.log2(q))
        return ans

    def crossEntropy(self, samples):
        """
        this measure computes the cross entropy with respect
        to some unknown probability distribution from which only samples
        are available. This measure is known to minimize the kl divergence.

        @param samples: numpy array
        """
        q = np.zeros(samples.shape[0])
        for i, sample in enumerate(samples):
            q[i] = max(1e-10, self.pdf(sample))

        # compute the cross entropy
        return -np.mean(np.log2(q))

    def l2error(self, dist, testSamplesUnit=None, testSamplesProb=None,
                n=1e4, dtype=SampleType.ACTIVEPROBABILISTIC):
        r"""
        mean squared error, defined as

        || p - p_n ||^2 = \int (p(x) - p_n(x))^2 * p(x) dx
                        ~ 1/n \sum_i (p(x_i) - p_n(x_i)^2
        for x_i drawn from p.
        @param dist: Dist
        @param testSamplesUnit: numpy array
        @param testSamplesProb: numpy array
        @param n: int, if no test samples are given, just select them
                  uniformly within he range of the distribution
        """
        if testSamplesUnit is None:
            # generate uniformly distributed samples
            testSamplesUnit = np.random.rand(dist.getDim() * n)
            testSamplesProb = np.array(testSamplesUnit)

            # scale them linearily to the right range
            for i, (a, b) in enumerate(dist.getBounds()):
                testSamplesUnit[:, i] = (b - a) * testSamplesUnit[:, i] + a

            for i, (a, b) in enumerate(self.getBounds()):
                testSamplesProb[:, i] = (b - a) * testSamplesProb[:, i] + a

            # change the sample type
            dtype = SampleType.ACTIVEUNIT

        if testSamplesProb is None:
            testSamplesProb = testSamplesUnit

        n = testSamplesUnit.shape[0]
        # compute the l2error
        err = 0.
        for i in xrange(n):
            erri = (self.pdf(testSamplesProb[i, :]) -
                    dist.pdf(testSamplesUnit[i, :])) ** 2

            if dtype == SampleType.ACTIVEUNIT:
                erri *= dist.pdf(testSamplesUnit[i, :])

            err += erri

        return err / n

    @classmethod
    def fromJson(cls, jsonObject):
        import pysgpp.extensions.datadriven.uq.dists as dists
        if 'uq.dists.J' in jsonObject['module']:
            return dists.J.fromJson(jsonObject)
        elif 'uq.dists.Uniform' in jsonObject['module']:
            return dists.Uniform.fromJson(jsonObject)
        elif 'uq.dists.TNormal' in jsonObject['module']:
            return dists.TNormal.fromJson(jsonObject)
        elif 'uq.dists.Normal' in jsonObject['module']:
            return dists.Normal.fromJson(jsonObject)
        elif 'uq.dists.Lognormal' in jsonObject['module']:
            return dists.Lognormal.fromJson(jsonObject)
        elif 'uq.dists.TLognormal' in jsonObject['module']:
            return dists.TLognormal.fromJson(jsonObject)
        elif 'uq.dists.Beta' in jsonObject['module']:
            return dists.Beta.fromJson(jsonObject)
        elif 'uq.dists.MultivariateNormal' in jsonObject['module']:
            return dists.MultivariateNormal.fromJson(jsonObject)
        elif 'uq.dists.SGDEdist' in jsonObject['module']:
            return dists.SGDEdist.fromJson(jsonObject)
        elif 'uq.dists.KDEDist' in jsonObject['module']:
            return dists.KDEDist.fromJson(jsonObject)
        elif 'uq.dists.Datadist' in jsonObject['module']:
            return dists.DataDist.fromJson(jsonObject)
        elif 'uq.dists.NatafDist' in jsonObject['module']:
            return dists.NatafDist.fromJson(jsonObject)
        else:
            raise TypeError('Unknown distribution "%s" => Please register it in fromJson function' % jsonObject['module'])
