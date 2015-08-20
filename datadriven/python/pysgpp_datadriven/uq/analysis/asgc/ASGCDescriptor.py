from ASGCAnalysisSpecification import ASGCAnalysisSpecification

from pysgpp.extensions.datadriven.uq.estimators import (
    MonteCarloStrategy,
    IntegralStrategy,
    CollocationPointsStrategy)

from pysgpp.extensions.datadriven.uq.refinement import RefinementDescriptor
# import types


class ASGCDescriptor(object):

    def __init__(self, builder):
        self._builder = builder
        self.__specification = ASGCAnalysisSpecification()

    def withParameters(self, params):
        self.__specification.setParameters(params)
        return self

    def withRefinement(self):
        return RefinementDescriptor(self)

    def withQoI(self, qoi):
        self.__specification.setQoI(qoi)
        return self

    def withMCEstimationStrategy(self, n=1000):
        self.__specification.setEstimationStrategy(MonteCarloStrategy(n))
        return self

    def withCPEstimationStrategy(self):
        self.__specification.setEstimationStrategy(CollocationPointsStrategy())
        return self

    def withIntegralEstimationStrategy(self):
        self.__specification.setEstimationStrategy(IntegralStrategy())
        return self

    def withNumberOfMoments(self, k):
        self.__specification.setNumberOfMoments(k)
        return self

    def withExpectationValue(self, E, t=-1):
        self.__specification.addReferenceValue('E', E, t)
        return self

    def withVariance(self, V, t=-1):
        self.__specification.addReferenceValue('V', V, t)
        return self

    def withAnalyticalSurface(self, f, t=-1):
        self.__specification.addReferenceValue('f', f, t)
        return self

    def withTimeStepsOfInterest(self, ts):
        # check for uniqueness of time steps
        if len(ts) != len(set(ts)):
            raise AttributeError('time steps of interest are not unique')
        self.__specification.setTimeStepsOfInterest(ts)
        return self

    def getSpecification(self):
        return self.__specification
