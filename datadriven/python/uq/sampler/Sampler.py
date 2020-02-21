# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

class Sampler(object):

    def __init__(self):
        self._params = None
        self._dim = None
        self._iteration = 0

    def setParameters(self, params):
        self._params = params
        self._dim = params.getStochasticDim()

    def nextSamples(self, *args, **kws):
        raise NotImplementedError()

    def learnData(self, *args, **kws):
        raise NotImplementedError()

    def hasMoreSamples(self):
        raise NotImplementedError()
