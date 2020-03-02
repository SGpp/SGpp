#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org
#
"""
@file    InverseCDFTransformation.py
@brief This class transforms a point in [0, 1]^d according to the
inverse of a given distrubution function

@version  0.1

"""

from pysgpp.extensions.datadriven.uq.transformation.Transformation import Transformation

from pysgpp.extensions.datadriven.uq import jsonLib as ju


class RosenblattTransformation(Transformation):
    """
    The inverse CDF transformation class
    """

    def __init__(self, dist):
        """
        Constructor
        """
        super(RosenblattTransformation, self).__init__()
        self._dist = dist

    def unitToProbabilistic(self, p):
        """
        Performs the PPF function of the given distribution on p which
        has to be in [0, 1]
        @param p: numeric value in [0, 1]
        @return: numeric value in the probabilistic space
        """
        return self._dist.ppf(p)

    def probabilisticToUnit(self, q):
        """
        Performs the CDF function of the given distribution on p
        @param q: numeric value in the probabilistic space
        @return: numeric value in [0, 1]
        """
        return self._dist.cdf(q)

    def vol(self):
        return 1.

    def getSize(self):
        return self._dist.getDim()

    def __str__(self):
        return "rosenblattTransformation: %s -> U(0, 1)" % str(dist)

    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        attrName = "_dist"
        attrValue = self._dist
        serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        return "{" + s + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        from pysgpp.extensions.datadriven.uq.dists.Dist import Dist

        key = '_dist'
        if key in jsonObject:
            dist = Dist.fromJson(jsonObject[key])

        return RosenblattTransformation(dist)
