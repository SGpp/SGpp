from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes


class ASGCSamplerSpecification(object):
    """
    Collection of parameters, which specify the ASGC process.
    An object of the class is aggregated by the ASGC object.
    """
    def __init__(self):
        self.__params = None
        self.__knowledgeTypes = [KnowledgeTypes.SIMPLE]
        self._qoi = '_'
        self.__timeStepsOfInterest = [0]
        self.__stopPolicy = None
        self.__testSet = None

        self.__learnWithTest = False
        self.__learnWithFolding = False

    def getTestSet(self):
        return self.__testSet

    def setTestSet(self, value):
        self.__testSet = value

    def getLearnWithFolding(self):
        return self.__learnWithFolding

    def setLearnWithFolding(self, value):
        self.__learnWithFolding = value

    def getParameters(self):
        return self.__params

    def setParameters(self, value):
        self.__params = value

    def getLearnWithTest(self):
        return self.__learnWithTest

    def setLearnWithTest(self, value):
        self.__learnWithTest = value

    def getStopPolicy(self):
        return self.__stopPolicy

    def setStopPolicy(self, value):
        self.__stopPolicy = value

    def getRefinement(self):
        return self._refinement

    def getKnowledgeTypes(self):
        return self.__knowledgeTypes

    def getQoI(self):
        return self._qoi

    def getTimeStepsOfInterest(self):
        return self.__timeStepsOfInterest

    def setRefinement(self, value):
        self._refinement = value

    def setKnowledgeTypes(self, value):
        self.__knowledgeTypes = value

    def setQoI(self, value):
        self._qoi = value

    def setTimeStepsOfInterest(self, value):
        self.__timeStepsOfInterest = value
