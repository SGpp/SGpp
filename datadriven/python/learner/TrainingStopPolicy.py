# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


import types

## The class implements heuristics for testing if the learn process should be finished before learner is overfitted.
#  
# The test is made by calling method <code>
# @link isTrainingComplete() isTrainingComplete(learner)@endlink</code> of the class, which 
# returns True if training process should be finished.
class TrainingStopPolicy(object):
    __adaptiveIterationLimit = None     #Maximal number of refinement iterations
    __epochsLimit = None                #Maximal number of iterations, during which accuracy can decreases
    __MSELimit = None                   #MSE on validation data, that have to be achieved
    __accuracyLimit = None              #accuracy on validation data, that have to be achieved
    __gridSizeLimit = None                   #Maximal grid size
    
    
    ##Contructor
    def __init__(self):
        self.__adaptiveIterationLimit = None
        self.__epochsLimit = None
        self.__MSELimit = None
        self.__accuracyLimit = None
        self.__gridSizeLimit = None
        self.__oldGridSize = 0


    ## Returns the maximal number of refinement iterations
    # @return: the maximal number of refinement iterations
    def getAdaptiveIterationLimit(self):
        return self.__adaptiveIterationLimit


    ## Returns the maximal number of iterations, during which accuracy can decreases
    # @return: the maximal number of iterations, during which accuracy can decreases
    def getEpochsLimit(self):
        return self.__epochsLimit


    ## Returns MSE on validation data, that have to be achieved
    # @return: MSE on validation data, that have to be achieved
    def getMSELimit(self):
        return self.__MSELimit


    ## Returns the accuracy on validation data, that have to be achieved
    # @return: accuracy on validation data, that have to be achieved
    def getAccuracyLimit(self):
        return self.__accuracyLimit


    ## Returns the maximal grid size
    # @return: maximal grid size
    def getGridSizeLimit(self):
        return self.__gridSizeLimit

    
    
    ## Checks if learning process have to be stopped
    # @param learner: Learner object 
    # @return: boolean value, true if learning has to stop, false otherwise
    def isTrainingComplete(self, learner):
        ans = self.hasLimitReached(learner) or not self.hasGridSizeChanged(learner)
        self.__oldGridSize = self.getGridSize(learner)
        return ans

    
    def hasGridSizeChanged(self, learner):
        return self.__oldGridSize != self.getGridSize(learner)
    
    
    def getGridSize(self, learner):
        try:
            grid_size = learner.grid_size
        except AttributeError:
            grid_size = learner.grid.getSize()
        return grid_size


    def hasLimitReached(self, learner):
        ans = (self.__adaptiveIterationLimit is None or
               self.__adaptiveIterationLimit <= learner.getCurrentIterationNumber()) \
               and (self.getGridSizeLimit() is None or
                    self.getGridSizeLimit() <= learner.grid.getSize()) \
               and (self.getMSELimit() is None or
                    self.getMSELimit() >= learner.trainingOverall[-1])
        return ans


    # # Setter for Maximal number of refinement iterations
    # @param limit: integer Maximal number of refinement iterations
    def setAdaptiveIterationLimit(self, limit):
        self.__adaptiveIterationLimit = limit


    ## Setter for epochs limit
    # @param limit: integer Maximal number of iterations, during which accuracy can decreases
    def setEpochsLimit(self, limit):
        self.__epochsLimit = limit

    
    ## Setter for MSE limit
    # @param limit: double minimal MSE on validation data, that have to be achieved
    def setMSELimit(self, limit):
        self.__MSELimit = limit

    
    ## Setter for accuracy limit
    # @param limit: double accuracy on validation data, that have to be achieved
    def setAccuracyLimit(self, limit):
        self.__accuracyLimit = limit
    
    
    ## Setter for maximal grid size
    # @param limit: integer maximal grid size  
    def setGridSizeLimit(self, limit):
        self.__gridSizeLimit = limit
        
        
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object. )
    def toString(self):
        serializationString = "'module' : '" + self.__module__ + "',\n"
        for attrName in dir(self):
            attrValue = self.__getattribute__(attrName)
            
            #lists, integers, floats, dictionaries can serialized with str()
            if type(attrValue) in [types.ListType, types.IntType, 
                             types.FloatType, types.DictType] and attrName.find("__") != 0: 
                serializationString += "'" + attrName + "'" + " : " + str(attrValue) + ",\n"
            # serialize strings with quotes    
            elif type(attrValue) == types.StringType and attrName.find("__") != 0:
                serializationString += "'" + attrName + "'" + " : '" + attrValue + "',\n"
                
        serializationString = "{" + serializationString.rstrip(",\n") + "}"
        return serializationString
    
    
    ## Restores the TrainingStopPolicy object from the json object with attributes.
    #
    # @param cls python keyword (do not specify)
    # @param jsonObject A json object.
    # @return The restored TrainingStopPolicy object.
    @classmethod
    def fromJson(cls, jsonObject):
        policy = TrainingStopPolicy()
        if jsonObject.has_key('_TrainingStopPolicy__adaptiveIterationLimit'):
            policy.setAdaptiveIterationLimit(jsonObject['_TrainingStopPolicy__adaptiveIterationLimit'])
        if jsonObject.has_key('_TrainingStopPolicy__epochsLimit'):
            policy.setEpochsLimit(jsonObject['_TrainingStopPolicy__epochsLimit'])
        if jsonObject.has_key('_TrainingStopPolicy__MSELimit'):
            policy.setMSELimit(jsonObject['_TrainingStopPolicy__MSELimit'])
        if jsonObject.has_key('_TrainingStopPolicy__accuracyLimit'):
            policy.setAccuracyLimit(jsonObject['_TrainingStopPolicy__accuracyLimit'])
        if jsonObject.has_key('_TrainingStopPolicy__gridSize'):
            policy.setGridSizeLimit(jsonObject['_TrainingStopPolicy__gridSize'])
        if jsonObject.has_key('_TrainingStopPolicy__oldGridSize'):
            policy.__oldGridSize = jsonObject['_TrainingStopPolicy__oldGridSize']
        return policy

            
        

