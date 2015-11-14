from pysgpp.extensions.datadriven.uq.estimators import (MonteCarloStrategy,
                               AnalyticEstimationStrategy,
                               CollocationPointsStrategy)

from ASGCAnalysis import ASGCAnalysis


class ASGCAnalysisBuilder(object):

    def __init__(self):
        """
        Default constructor
        """
        self.__learner = None
        self.__estimationStrategy = None

    def withLearner(self, learner):
        self.__learner = learner
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
        if self.__learner is None:
            raise AttributeError('no learner specified')

        return ASGCAnalysis(self.__learner, strategy=self.__estimationStrategy)
