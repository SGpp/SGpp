from pysgpp.extensions.datadriven.uq.dists import Dist


# # Collection of parameters, which specify the ASGC process.
# An object of the class is aggregated by the ASGC object.
class ASGCAnalysisSpecification(object):

    def __init__(self):
        self.__params = None
        self.__distribution = None
        self.__estimator = None
        self._refinement = None
        self.__estimationStrategy = None
        self.__k = None

        self._qoi = '_'

        self.__reachesSteadyState = False
        self.__timeStepsOfInterest = [0]

        self.__referenceValues = {}

    def setTimeStepsOfInterest(self, ts):
        self.__timeStepsOfInterest = ts

    def getTimeStepsOfInterest(self):
        return self.__timeStepsOfInterest

    def setRefinement(self, refinement):
        self._refinement = refinement

    def getRefinement(self):
        return self._refinement

    def setEstimator(self, estimator):
        self.__estimator = estimator

    def getEstimator(self):
        return self.__estimator

    def setEstimationStrategy(self, strategy):
        self.__estimationStrategy = strategy

    def getEstimationStrategy(self):
        return self.__estimationStrategy

    def setParameters(self, params):
        self.__params = params
        self.__distribution = params.getIndependentJointDistribution()

    def getParameters(self):
        return self.__params

    def setDistribution(self, distribution):
        if issubclass(distribution.__class__, Dist):
            self.__distribution = distribution

    def getDistribution(self):
        return self.__distribution

    def getDim(self):
        return self.__params.getDim()

    def setQoI(self, qoi):
        self._qoi = qoi

    def getQoI(self):
        return self._qoi

    def setNumberOfMoments(self, k):
        self.__k = k

    def getNumberOfMoments(self):
        return self.__k

    def addReferenceValue(self, key, value, t=-1):
        if t not in self.__referenceValues:
            self.__referenceValues[t] = {}

        self.__referenceValues[t][key] = value

    def getReferenceValues(self):
        return self.__referenceValues
