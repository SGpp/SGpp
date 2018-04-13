#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    Linear.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 13:38:59 2013

@brief The class performs a linear transformation from the unit
hypercube to an arbitrary hypercube

@version  0.1
"""
from Transformation import Transformation

import pysgpp.extensions.datadriven.uq.jsonLib as ju


class LinearTransformation(Transformation):

    def __init__(self, lower, upper):
        self.__a, self.__b = lower, upper

    def unitToProbabilistic(self, x):
        """
        Performs a linear transformation of x in [0, 1] to [a, b]
        @param x: float value
        @return: transformed value
        """
        if self.__a == 0 and self.__b == 1:
            return x
        elif self.__a == self.__b:
            return self.__a
        else:
            return 1. * x * (self.__b - self.__a) + self.__a

    def probabilisticToUnit(self, x):
        """
        Performs a linear transformation of x in [a, b] to [0, 1]
        @param x: float value
        @return: transformed value
        """
        if self.__a == 0 and self.__b == 1:
            return x
        elif self.__a == self.__b:
            return self.__a
        else:
            return (x - self.__a) / (self.__b - self.__a)

    def vol(self):
        return self.__b - self.__a

    def getBounds(self):
        return [self.__a, self.__b]

    def getSize(self):
        return 1

    def __str__(self):
        return "LinearTransformation: [%g, %g] -> [%g, %g]" % (self.__a,
                                                               self.__b,
                                                               0, 1)

    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName in ("_LinearTransformation__a", "_LinearTransformation__b"):
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        return "{" + s + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        key = '_LinearTransformation__a'
        if key in jsonObject:
            a = float(jsonObject[key])

        key = '_LinearTransformation__b'
        if key in jsonObject:
            b = float(jsonObject[key])

        return LinearTransformation(a, b)
