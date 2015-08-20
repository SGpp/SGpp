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

    def __init__(self, samplePath):
        """
        Constructor. There are some restrictions to the samples:
        As they represent the underlying probability, they have to be stored
        in the probabilistic space and within each other they are uniformly
        distributed.
        @param samplePath: path to the file where the samples are stored.
        """
        # check if the data exists
        if samplePath is None or not os.path.exists(samplePath):
            raise AttributeError('the file "%s" does not exist' % samplePath)

        # read data and store it in a set of samples
        self.samples = readDataTrivial(samplePath, ' ', False)['data'].array()
        self.__sampleToIndex = {}
        for sample in self.samples:
            self.__sampleToIndex[tuple(sample)] = len(self.samples)

        # store number of samples and their dimensionality
        self.__n, self.__dim = self.samples.shape

        # find the bounds of the data
        mins = np.apply_along_axis(np.min, 1, self.samples)
        maxs = np.apply_along_axis(np.max, 1, self.samples)
        self.__bounds = np.vstack((mins, maxs)).T

    def pdf(self, p, *args, **kws):
        if tuple(p) in self.__sampleToIndex:
            return 1. / len(self.samples)
        else:
            return 0.

    def cdf(self, p, *args, **kws):
        ans = np.ndarray(self.__dim, dtype='int')
        for i, pi in enumerate(p):
            ans[i] = np.sum([1 for qi in self.samples[:, i] if pi <= qi])
        return ans

    def ppf(self, p, *args, **kws):
        raise NotImplementedError()

    def rvs(self, n=1):
        ixs = np.random.randint(0, len(self.samples) - 1, n)
        return self.samples[ixs, :]

    def mean(self):
        return np.mean(self.samples, axis=0, dtype='float32')

    def var(self):
        return np.var(self.samples, ddof=1, axis=0, dtype='float32')

    def std(self):
        return np.std(self.samples, ddof=1, axis=0, dtype='float32')

    def getBounds(self):
        return self.__bounds

    def getDim(self):
        return self.__dim

    @classmethod
    def fromJson(cls, jsonObject):
        return Dist.fromJson(cls, jsonObject)
