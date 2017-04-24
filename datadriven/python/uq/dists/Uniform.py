#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    uniform.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 14:26:01 2013

@brief   uniform distribution

@version  0.1

"""


from Dist import Dist
from scipy.stats import uniform
import numpy as np

import pysgpp.extensions.datadriven.uq.jsonLib as ju


class Uniform(Dist):
    """
    Represents a uniform distribution
    """

    def __init__(self, a, b):
        """
        Constructor
        @param a: lower interval threshold
        @param b: upper interval threshold
        """
        super(Uniform, self).__init__()
        self.__a = float(a)
        self.__b = float(b)

        if a >= b:
            raise AttributeError('lower bound of the interval is larger then \
                                  the higher one')

        self._dist = uniform(loc=a, scale=b-a)

    def pdf(self, x):
        return self._dist.pdf(x)

    def cdf(self, x):
        return self._dist.cdf(x)

    def ppf(self, x):
        return self._dist.ppf(x)

    def rvs(self, n=1):
        return self._dist.rvs(n)

    def mean(self):
        return self._dist.mean()

    def var(self):
        return self._dist.var()

    def std(self):
        return self._dist.std()

    def getBounds(self):
        return np.array([self.__a, self.__b], dtype="float")

    def getDim(self):
        return 1

    def __str__(self):
        return "U(%g, %g)" % (self.__a, self.__b)

    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName in ("_Uniform__a", "_Uniform__b"):
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        return "{" + s + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the Uniform object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored Uniform object
        """
        key = '_Uniform__a'
        if key in jsonObject:
            a = jsonObject[key]

        key = '_Uniform__b'
        if key in jsonObject:
            b = jsonObject[key]

        return Uniform(a, b)
