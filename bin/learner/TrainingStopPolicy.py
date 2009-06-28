#!/usr/bin/env python
# -*- coding: utf-8 -*-


class TrainingStopPolicy(object):
    """ generated source for TrainingStopPolicy

    """
    __adaptiveIterationLimit = None
    __epochsLimit = None
    __MSELimit = None
    __accuracyLimit = None
    __gridSize = None
    
    def __init__(self):
        self.__adaptiveIterationLimit = None
        self.__epochsLimit = None
        self.__MSELimit = None
        self.__accuracyLimit = None
        self.__gridSize = None

    #@todo: Make it more advanced
    def isTrainingComplete(self, learner):
        if self.__adaptiveIterationLimit != None and self.__adaptiveIterationLimit < learner.getCurrentIterationNumber():
            return True
        return False

    def setAdaptiveIterationLimit(self, limit):
        self.__adaptiveIterationLimit = limit

    def setEpochsLimit(self, limit):
        self.__epochsLimit = limit

    def setMSELimit(self, limit):
        self.__MSELimit = limit

    def setAccuracyLimit(self, limit):
        self.__accuracyLimit = limit
        
    def setGridSizeLimit(self, limit):
        self.__gridSize = limit
        


