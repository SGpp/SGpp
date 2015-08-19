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

from scipy.stats import multivariate_normal
from Dist import Dist
import numpy as np
import pysgpp_datadriven.uq.jsonLib as ju


class MultivariateNormal(Dist):
    """
    Represents a multivariate normal distribution
    """

    def __init__(self, mu, cov, a, b):
        """
        Constructor
        @param mu: mean of 1d gaussians
        @param cov: covariance matrix
        @param a: lower boundary
        @param b: upper boundary
        """
        super(MultivariateNormal, self).__init__()
        self.__mu = mu
        self.__cov = cov
        self.__a, self.__b = a, b
        self.__dim = len(mu)

        # standard multivariate normal
        self._dist = multivariate_normal(mu, cov)

    def pdf(self, x):
        return self._dist.pdf(x)

    def rvs(self, n=1):
        return self._dist.rvs(n)

    def getBounds(self):
        ans = np.zeros([self.__dim, 2], dtype="float")
        ans[:, 0] = self.__a
        ans[:, 1] = self.__b
        return ans

    def getDim(self):
        return self.__dim

    def __str__(self):
        return "mN^(%i, %i, %g, %g)" % (self.__dim, self.__dim,
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
        key = '_MultivariateNormal__mu'
        if key in jsonObject:
            mu = jsonObject[key]

        key = '_MultivariateNormal__cov'
        if key in jsonObject:
            cov = jsonObject[key]

        key = '_MultivariateNormal__a'
        if key in jsonObject:
            a = jsonObject[key]

        key = '_MultivariateNormal__b'
        if key in jsonObject:
            b = jsonObject[key]

        return MultivariateNormal(mu, cov, a, b)
