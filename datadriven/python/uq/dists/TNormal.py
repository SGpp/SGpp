#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    tnormal.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 14:22:17 2013

@brief   truncated normal distribution by some confidence value alpha

@version  0.1

"""

from scipy.stats import truncnorm, norm
from Dist import Dist

import pysgpp.extensions.datadriven.uq.jsonLib as ju

import numpy as np


class TNormal(Dist):
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
        super(TNormal, self).__init__()
        self.__mu = float(mu)
        self.__sigma = float(sigma)
        self.__a, self.__b = a, b

        # truncated standard normal
        a, b = self.__trans(a), self.__trans(b)
        self._dist = truncnorm(a, b)

    @classmethod
    def by_range(cls, *args, **kws):
        """
        Constructor given a interval
        """
        cls(*args, **kws)

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

    def __trans(self, x):
        return (x - self.__mu) / self.__sigma

    def __inv_trans(self, y):
        return y * self.__sigma + self.__mu

    def pdf(self, x):
        return self._dist.pdf(self.__trans(x)) / self.__sigma

    def cdf(self, x):
        return self._dist.cdf(self.__trans(x))

    def ppf(self, x):
        return self.__inv_trans(self._dist.ppf(x))

    def mean(self):
        return self.__mu

    def var(self):
        return self.__sigma * self.__sigma

    def std(self):
        return self.__sigma

    def rvs(self, n=1):
        return self.__inv_trans(self._dist.rvs(n))

    def getBounds(self):
        return np.array([self.__a, self.__b], dtype='float')

    def getDim(self):
        return 1

    def __str__(self):
        return "TN(%g, %g, %g, %g)" % (self.__mu, self.__sigma,
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
        key = '_TNormal__my'
        if key in jsonObject:
            my = jsonObject[key]

        key = '_TNormal__sigma'
        if key in jsonObject:
            sigma = jsonObject[key]

        key = '_TNormal__xlim'
        if key in jsonObject:
            xlim = jsonObject[key]

        return TNormal(my, sigma, xlim[0], xlim[1])
