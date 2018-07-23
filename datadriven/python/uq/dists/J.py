#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    j.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 14:14:52 2013

@brief   Joint multivariate distribution of independent distributions

@version  0.1

"""
from Dist import Dist
from SGDEdist import SGDEdist
import numpy as np
from pysgpp.extensions.datadriven.uq.operations import discretizeFunction


class J(Dist):
    """
    Models the multivariate distribution of independent variables
    """

    def __init__(self, dists):
        super(J, self).__init__()
        if all([issubclass(dist.__class__, Dist) for dist in dists]):
            self.__dists = np.array(dists)
            self.__ixs = []
            self.__dim = 0
            self.__n = len(self.__dists)
            i = 0
            for dist in self.__dists:
                di = dist.getDim()
                self.__dim += di
                self.__ixs += [[i + xi for xi in range(di)]]
                i += di
            # print self.__dim, self.__ixs
            self.sgdeDist = None
            self.error = None
        else:
            raise AttributeError('Not valid distributions')

    def pdf(self, p, marginal=False):
        ans = np.ndarray(self.__n, dtype='float')
        for i, ix in enumerate(self.__ixs):
            if len(ix) == 1:
                x = p[ix[0]]
            else:
                x = np.array([p[j] for j in ix], dtype='float')
            ans[i] = self.__dists[i].pdf(x)

        if marginal:
            return ans
        else:
            return np.prod(ans, dtype='float')

    def cdf(self, p):
        ans = np.ndarray(self.__n, dtype='float')
        for i, ix in enumerate(self.__ixs):
            if len(ix) == 1:
                x = p[ix[0]]
            else:
                x = np.array([p[j] for j in ix], dtype='float')
            ans[i] = self.__dists[i].cdf(x)
        return ans

    def ppf(self, p):
        ans = np.ndarray(self.__n, dtype='float')
        for i, ix in enumerate(self.__ixs):
            if len(ix) == 1:
                x = p[ix[0]]
            else:
                x = np.array([p[j] for j in ix], dtype='float')
            ans[i] = self.__dists[i].ppf(x)
        return ans

    def mean(self):
        m = np.array([d.mean() for d in self.__dists])
        return np.prod(m)

    def var(self):
        a = np.array([d.var() for d in self.__dists])
        b = np.array([d.mean() for d in self.__dists])
        return np.prod(a + b ** 2) - np.prod(b ** 2)

    def std(self):
        return np.sqrt(self.var())

    def rvs(self, n=1):
        ans = np.ndarray([self.__dim, n], dtype='float')
        for i, ix in enumerate(self.__ixs):
            r = self.__dists[i].rvs(n)
            # if r is a matrix
            if len(ix) > 1:
                for j in xrange(len(ix)):
                    ans[ix[j], :] = r[:, j]
            # ... if it is just a vector
            else:
                ans[ix[0], :] = r

        return ans.T

    def getBounds(self):
        ans = np.ndarray([self.__dim, 2], dtype='float')
        for i, ix in enumerate(self.__ixs):
            ans[ix, :] = self.__dists[i].getBounds()
        return ans

    def getDistributions(self):
        return [d for d in self.__dists]

    def getTupleIndices(self):
        return self.__ixs

    def __getitem__(self, index):
        return self.__dists[index]

    def getDim(self):
        return self.__dim

    def discretize(self, level=5, *args, **kws):
        """
        discretize the pdf of the current distribution
        using a sparse grid interpolant
        """
        if self.sgdeDist is None:
            def f(p):
                return float(self.pdf(p))

            # discretize pdf
            grid, alpha, self.error = \
                discretizeFunction(f, self.getBounds(),
                                   level=level, *args, **kws)
            # generate a SGDE function
            self.sgdeDist = SGDEdist(grid, alpha)

        return self.sgdeDist, self.error

#
#     def discretize(self, *args, **kws):
#         ans = [None] * self.__n
#         for i, dist in enumerate(self.__dists):
#             ans[i] = dist.discretize(*args, **kws)
#         if len(ans) > 1:
#             grid, _, error = ans[0]
#             err = error[-1, 1]
#             for i in xrange(1, len(ans)):
#                 # join the grids
#                 grid2, _, error2 = ans[i]
#                 grid = join(grid, grid2)
#
#                 # accumulate the error
#                 err += error2[-1, 1]
#
#         # hierarchize the result
#         gs = grid.getStorage()
#         p = DataVector(gs.getDimension())
#         nodalValues = DataVector(gs.size())
#         for i in xrange(gs.size()):
#             gs.getPoint(i).getStandardCoordinates(p)
#             nodalValues[i] = self.pdf(p.array())
#
#         alpha = hierarchize(grid, nodalValues)
#
#         return grid, alpha, err

    def __str__(self):
        return "J(%s)" % str([str(dist) for dist in self.__dists])

    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        # serialize dists
        attrName = "_J__dists"
        attrValue = self.__getattribute__(attrName)
        x = [dist.toJson() for dist in attrValue]
        x = ['"' + str(i) + '": ' + str(xi) for i, xi in enumerate(x)]
        serializationString += '"' + attrName + '": {' + ', '.join(x) + '}'

        return "{" + serializationString + "} \n"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the J object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored J object
        """
        key = '_J__dists'
        if key in jsonObject:
            vals = jsonObject[key]
            dists = [Dist.fromJson(vals[key]) for key in sorted(vals.keys())]
        else:
            raise AttributeError("J: fromJson - the mandatory keyword '%s' does not exist" % key)

        return J(dists)
