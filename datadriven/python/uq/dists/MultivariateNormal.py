#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
from pysgpp.extensions.datadriven.uq.dists.Normal import Normal
"""
@file    tnormal.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 14:22:17 2013

@brief   truncated normal distribution by some confidence value alpha

@version  0.1

"""

from J import J
from Dist import Dist
import numpy as np
import pysgpp.extensions.datadriven.uq.jsonLib as ju

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
        dists = [None] * self.__dim
        for idim in xrange(self.__dim):
            dists[idim] = Normal.by_alpha(0, 1, 0.001)
        self.dist = J(dists)

        self.cov_inv = np.linalg.inv(cov)
        self.norm = 1. / np.sqrt((2 * np.pi) ** self.__dim * np.linalg.det(cov))

        # do cholesky decomposition for nataf transformation
        self.corr = np.ndarray(cov.shape)
        for i in xrange(cov.shape[0]):
            for j in xrange(cov.shape[1]):
                self.corr[i, j] = cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])
        self.L = np.linalg.cholesky(self.corr)

    def pdf(self, x):
        if all([self.__a <= xi <= self.__b for xi in x]):
            z = x - self.__mu
            return self.norm * np.exp(-0.5 * np.dot(z, np.dot(self.cov_inv, z)))
        else:
            return 0.0

    def rvs(self, n=1):
        # do a nataf transformation
        isample = 0
        ans = np.ndarray((n, self.__dim))
        while (isample < n):
            sample = np.dot(self.L, self.dist.rvs(1).reshape((self.__dim, 1))).T[0]
            for idim in xrange(self.__dim):
                sample[idim] = self.__mu[idim] + sample[idim] * np.sqrt(self.__cov[idim, idim])
            if self.__a < sample.min() and sample.max() < self.__b:
                ans[isample, :] = sample
                isample += 1
        return ans

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
