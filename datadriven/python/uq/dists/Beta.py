#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org
#
"""
@file    beta.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 14:10:32 2013

@brief   Represents the beta distribution

@version  0.1

"""

from Dist import Dist
import pysgpp.extensions.datadriven.uq.jsonLib as ju
import scipy.stats as stats


class Beta(Dist):
    """
    The Beta distribution class
    """

    def __init__(self, p, q, accLevel=0., width=1.):
        """
        Constructor
        @param p: first design parameter
        @param q: second design parameter
        @param accLevel: lower boundary
        @param width: width of the interval
        """
        super(Beta, self).__init__()
        self.__p = p
        self.__q = q
        self._dist = stats.beta(p, q, accLevel, width)

        self.__l = accLevel
        self.__width = width

    def pdf(self, x):
        return self._dist.pdf(x)

    def cdf(self, x):
        return self._dist.cdf(x)

    def ppf(self, x):
        return self._dist.ppf(x)

    def mean(self):
        return self._dist.mean()

    def var(self):
        return self._dist.var()

    def std(self):
        return self._dist.std()

    def rvs(self, n=1):
        return self._dist.rvs(n)

    def getBounds(self):
        return [self.__l, self.__l + self.__width]

    def getDim(self):
        return 1

    def alpha(self):
        return self.__p

    def beta(self):
        return self.__q

    def __str__(self):
        return "B(%i, %i, %g, %g)" % (self.__p, self.__q,
                                      self.__l, self.__width)

    def toJson(self):
        """
        Returns a string that represents the object

        Arguments:

        Return A string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName in ["_Beta__p", "_Beta__q", "_Beta__l", "_Beta__width"]:
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        return "{" + s + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the Beta object from the json object with its
        attributes.

        Arguments:
        jsonObject -- json object

        Return the restored UQSetting object
        """
        # restore surplusses
        key = '_Beta__p'
        if key in jsonObject:
            p = int(jsonObject[key])

        key = '_Beta__q'
        if key in jsonObject:
            q = int(jsonObject[key])

        key = '_Beta__l'
        if key in jsonObject:
            l = float(jsonObject[key])

        key = '_Beta__width'
        if key in jsonObject:
            width = float(jsonObject[key])

        return Beta(p, q, l, width)
