#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    Corr.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 14:14:52 2013

@brief   Correlated distribution

@version  0.1

"""

from Dist import Dist

import numpy as np


class Corr(Dist):
    """
    Models the multivariate distribution of correlated
    """

    def __init__(self, dists):
        super(Corr, self).__init__()
        if all([issubclass(dist.__class__, Dist) for dist in dists]):
            self.__dists = dists
            self.__dim = len(dists)
        else:
            raise AttributeError('Not valid distributions')

    def pdf(self, x):
        """
        When x is correlated to y, but y not to x then the joint
        probability density is defined as

        p(x, y) = p(x|y) * p(y)

        The correlation is described as a tuple (x, (y, )).
        """
        if len(x) != self.__dim:
            raise AttributeError('Expected tuple of length %i but got of \
                                  length %i' % (self.__dim, len(x)))

        p1 = self.__dists[0].pdf(x[0])
        p2 = self.__dists[1].pdf(x)

        return p1 * p2

    def rvs(self, n=1):
        ans = np.array([[0, 0]] * n, dtype='float')
        for i in xrange(n):
            x1 = self.__dists[0].rvs(1)[0]
            x2 = self.__dists[1].rvs(x1, 1)[0]
            ans[i, :] = [x1, x2]
        return ans

    def getBounds(self):
        ans = []
        for dist in self.__dists:
            bounds = dist.getBounds()
            if isinstance(bounds[0], list):
                for bound in bounds:
                    ans += [bound]
            else:
                ans += [bounds]
        return ans

    def getDistributions(self):
        return [d for d in self.__dists]

    def getDim(self):
        return self.__dim

    def __str__(self):
        return str([str(dist) for dist in self.__dists])

    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        # serialize dists
        attrName = "_Corr__dists"
        attrValue = self.__getattribute__(attrName)
        x = [dist.toJson() for dist in attrValue]
        x = ['"' + str(i) + '": ' + str(xi) for i, xi in enumerate(x)]
        serializationString += '"' + attrName + '": {' + ', '.join(x) + '}'

        return "{" + serializationString + "} \n"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the Corr object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored Corr object
        """
        key = '_Corr__dists'
        if key in jsonObject:
            vals = jsonObject[key]
            dists = [Dist.fromJson(vals[key]) for key in sorted(vals.keys())]

        return Corr(dists)
