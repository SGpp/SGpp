from bin.uq.estimators import (MonteCarloStrategy,
                               AnalyticEstimationStrategy,
                               CollocationPointsStrategy)

from ASGCAnalysis import ASGCAnalysis
from ASGCKnowledge import ASGCKnowledge
from ASGCKnowledgeFormatter import ASGCKnowledgeFormatter


class ASGCAnalysisBuilder(object):

    def __init__(self):
        """
        Default constructor
        """
        self.__learner = None
        self.__qoi = None
        self.__strategy = MonteCarloStrategy()

    def withLearner(self, learner):
        self.__learner = learner
        return self

    def withMonteCarloEstimationStrategy(self, *args, **kws):
        self.__strategy = MonteCarloStrategy(*args, **kws)
        return self

    def withAnalyticEstimationStrategy(self, *args, **kws):
        self.__strategy = AnalyticEstimationStrategy()
        return self

    def withCollocationPointsStrategy(self):
        self.__strategy = CollocationPointsStrategy()
        return self

    def andGetResult(self):
        """
        Returns the adaptive sparse grid collocation object that is
        currently being constructed
        """
        if self.__learner is None:
            raise AttributeError('no learner specified')

        return ASGCAnalysis(self.__learner, strategy=self.__strategy)
