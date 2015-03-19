#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org
#
"""
@file    InverseCDFTransformation.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 13:34:25 2013

@brief This class transforms a point in [0, 1]^d according to the
inverse of a given distrubution function

@version  0.1

"""

from Transformation import Transformation


class InverseCDFTransformation(Transformation):
    """
    The inverse CDF transformation class
    """

    def __init__(self, dist):
        """
        Constructor
        """
        super(InverseCDFTransformation, self).__init__()
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

    @classmethod
    def fromJson(cls, jsonObject):
        from bin.uq.dists.Dist import Dist

        key = '_InverseCDFTransformation__dist'
        if key in jsonObject:
            dist = Dist.fromJson(jsonObject[key])

        return InverseCDFTransformation(dist)
