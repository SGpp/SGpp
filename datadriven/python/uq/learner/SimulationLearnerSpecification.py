# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.learner.TrainingSpecification import TrainingSpecification
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes


class SimulationLearnerSpecification(TrainingSpecification):
    """
    Collection of parameters, which specify the ASGC process.
    An object of the class is aggregated by the ASGC object.
    """
    def __init__(self):
        super(self.__class__, self).__init__()
        # simulation based parameters
        self.__params = None
        self._qoi = '_'
        self.__knowledgeTypes = [KnowledgeTypes.SIMPLE]
        self.__timeStepsOfInterest = [0]

        # adaptive refinement
        self._refinement = None

    def getParameters(self):
        return self.__params

    def setParameters(self, value):
        self.__params = value

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

    def getDim(self):
        return self.__params.getDim()
