#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    Data.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Tue Jul 08 14:04:16 2014

@brief   discrete distribution based on data

@version  0.1
"""

from pysgpp.extensions.datadriven.tools import readDataTrivial
import os

from Dist import Dist
import numpy as np


class DataDist(Dist):
    """
    Models a discrete distribution given by data
    """

    def __init__(self, samples):
        """
        Constructor. There are some restrictions to the samples:
        As they represent the underlying probability, they have to be drawn
        iid.
        @param samples: numpy array (num_samples x num_dims)
        """
        # read data and store it in a set of samples
        self.samples = samples
        self.__sampleToIndex = {}
        for sample in self.samples:
            x = tuple(sample)
            if x in self.__sampleToIndex:
                self.__sampleToIndex[x] += 1
            else:
                self.__sampleToIndex[x] = 1

        # store number of samples and their dimensionality
        self.__n, self.__dim = self.samples.shape

        # find the bounds of the data
        mins = np.min(self.samples, axis=0)
        maxs = np.max(self.samples, axis=0)
        self.__bounds = np.vstack((mins, maxs)).T

    def pdf(self, p, *args, **kws):
        x = tuple(p)
        if x in self.__sampleToIndex:
            return float(self.__sampleToIndex[x]) / self.__n
        else:
            return 0.

    def cdf(self, p, *args, **kws):
        ans = np.ndarray(self.__dim, dtype='int')
        for i, pi in enumerate(p):
            ans[i] = np.sum([1. for qi in self.samples[:, i] if pi <= qi])
        return ans

    def ppf(self, p, *args, **kws):
        raise NotImplementedError()

    def rvs(self, n=1):
        ixs = np.random.randint(0, len(self.samples) - 1, n)
        return self.samples[ixs, :]

    def mean(self):
        return np.mean(self.samples, axis=0)

    def var(self):
        return np.var(self.samples, ddof=1, axis=0)

    def std(self):
        return np.sqrt(self.var())

    def getBounds(self):
        return self.__bounds

    def getDim(self):
        return self.__dim

    def __str__(self):
        return "Data((%i x %i), %s)" % (self.__n, self.__dim, self.__bounds)

    @classmethod
    def fromJson(cls, jsonObject):
        return Dist.fromJson(cls, jsonObject)
