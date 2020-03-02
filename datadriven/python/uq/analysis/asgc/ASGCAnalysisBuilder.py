# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.estimators import (MonteCarloStrategy,
                                                        AnalyticEstimationStrategy,
                                                        CollocationPointsStrategy)

from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCAnalysis import ASGCAnalysis


class ASGCAnalysisBuilder(object):

    def __init__(self):
        """
        Default constructor
        """
        self.__uqManager = None
        self.__estimationStrategy = None

    def withUQManager(self, uqManager):
        self.__uqManager = uqManager
        return self

    def withMonteCarloEstimationStrategy(self, *args, **kws):
        self.__estimationStrategy = MonteCarloStrategy(*args, **kws)
        return self

    def withAnalyticEstimationStrategy(self, *args, **kws):
        self.__estimationStrategy = AnalyticEstimationStrategy()
        return self

    def withCollocationPointsStrategy(self):
        self.__estimationStrategy = CollocationPointsStrategy()
        return self

    def andGetResult(self):
        """
        Returns the adaptive sparse grid collocation object that is
        currently being constructed
        """
        if self.__uqManager is None:
            raise AttributeError('no learner specified')

        return ASGCAnalysis(self.__uqManager, strategy=self.__estimationStrategy)
