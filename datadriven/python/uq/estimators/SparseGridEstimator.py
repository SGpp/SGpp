# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.estimators.MonteCarloStrategy import MonteCarloStrategy
from pysgpp.extensions.datadriven.uq.estimators.Estimator import Estimator


class SparseGridEstimator(Estimator):

    def __init__(self, params, strategy=None):
        """
        Constructor
        @param params: parameter set
        @param strategy: estimation strategy
        """
        self.__params = params.activeParams()
        self.__T = self.__params.getJointTransformation()
        self.__vol = self.__T.vol()
        self.__U = self.__params.getIndependentJointDistribution()
        if strategy is None:
            self.__estimationStrategy = MonteCarloStrategy()
        else:
            self.__estimationStrategy = strategy

    def mean(self, grid, alpha, *args, **kws):
        """
        Estimate moments of the given sparse grid function
        @param grid: Grid
        @param alpha: coefficients vector
        @return: tuple(moment, error)
        """
        return self.__estimationStrategy.mean(self.__vol, grid, alpha,
                                              self.__U, self.__T,
                                              *args, **kws)

    def var(self, grid, alpha, *args, **kws):
        """
        Estimate moments of the given sparse grid function
        @param grid: Grid
        @param alpha: coefficients vector
        @return: tuple(moment, error)
        """
        return self.__estimationStrategy.mean(self.__vol, grid, alpha,
                                              self.__U, self.__T,
                                              *args, **kws)

    def setStrategy(self, strategy):
        self.__estimationStrategy = strategy

    def getStrategy(self):
        return self.__estimationStrategy
