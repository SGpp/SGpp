#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
from pysgpp.extensions.datadriven.uq.transformation import LinearTransformation
"""
@file    tnormal.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 14:22:17 2013

@brief   truncated normal distribution by some confidence value alpha

@version  0.1

"""

from scipy.stats import norm
from Dist import Dist

import pysgpp.extensions.datadriven.uq.jsonLib as ju

import numpy as np


class Normal(Dist):
    """
    Represents a truncated normal distribution

    See: http://en.wikipedia.org/wiki/Truncated_normal_distribution
    """

    def __init__(self, mu, sigma, a, b):
        """
        Constructor
        @param mu: expectation value
        @param sigma: standard deviation
        @param a: lower boundary
        @param b: upper boundary
        """
        super(Normal, self).__init__()
        self.__mu = float(mu)
        self.__sigma = float(sigma)
        self.__a, self.__b = a, b

        # standard normal
        self._dist = norm(loc=mu, scale=sigma)

        self.__linearTrans = LinearTransformation(self._dist.cdf(self.__a),
                                                  self._dist.cdf(self.__b))

    @classmethod
    def by_range(cls, *args, **kws):
        """
        Constructor given a interval
        """
        return cls(*args, **kws)

    @classmethod
    def by_alpha(cls, mu, sigma, alpha):
        """
        Constructor given a confidence value
        @param mu: expectation value
        @param sigma: standard deviation
        @param alpha: confidence value
        """
        U = norm(loc=mu, scale=sigma)
        a = U.ppf(alpha / 2.)
        b = U.ppf(1 - alpha / 2.)

        return cls(mu, sigma, a, b)

    def pdf(self, x):
        if self.__a <= x <= self.__b:
            return self._dist.pdf(x)
        else:
            return 0.0;

    def cdf(self, x):
        if self.__a <= x <= self.__b:
            x_unit = self._dist.cdf(x)
            return self.__linearTrans.probabilisticToUnit(x_unit)
        else:
            raise AttributeError("normal: cdf - x out of range [%g, %g]" % (self.__a, self.__b))

    def ppf(self, x):
        if 0.0 <= x <= 1.0:
            x_unit = self.__linearTrans.unitToProbabilistic(x)
            return self._dist.ppf(x_unit)
        else:
            raise AttributeError("normal: ppf - x out of range [%g, %g]" % (self.__a, self.__b))

    def mean(self):
        return self.__mu

    def var(self):
        return self.__sigma * self.__sigma

    def std(self):
        return self.__sigma

    def rvs(self, n=1):
        samples = np.zeros(n)
        i = 0
        while i < n:
            newSamples = self._dist.rvs(n - i)
            # check range
            for sample in newSamples:
                if self.__a <= sample <= self.__b:
                    samples[i] = sample
                    i += 1
        return samples

    def getBounds(self):
        return np.array([self.__a, self.__b], dtype='float')

    def getDim(self):
        return 1

    def __str__(self):
        return "N(%g, %g, %g, %g)" % (self.__mu, self.__sigma,
                                      self.__a, self.__b)

    def toJson(self):
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName in dir(self):
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        return "{" + s + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the TNormal object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored TNormal object
        """
        key = '_Normal__mu'
        if key in jsonObject:
            mu = jsonObject[key]

        key = '_Normal__sigma'
        if key in jsonObject:
            sigma = jsonObject[key]

        key = '_Normal__a'
        if key in jsonObject:
            a = jsonObject[key]

        key = '_Normal__b'
        if key in jsonObject:
            b = jsonObject[key]

        return Normal(mu, sigma, a, b)
