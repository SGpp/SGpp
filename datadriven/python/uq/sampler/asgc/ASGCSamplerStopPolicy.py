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
class ASGCSamplerStopPolicy(object):
    ##Contructor
    def __init__(self):
        self.__adaptiveIterationLimit = None  # Maximal number of refinement iterations
        self.__accuracyLimit = None  # accuracy on validation data, that have to be achieved
        self.__gridSizeLimit = None  # Maximal grid size
        self.__oldGridSize = 0


    ## Returns the maximal number of refinement iterations
    # @return: the maximal number of refinement iterations
    def getAdaptiveIterationLimit(self):
        return self.__adaptiveIterationLimit


    ## Returns the accuracy on validation data, that have to be achieved
    # @return: accuracy on validation data, that have to be achieved
    def getAccuracyLimit(self):
        return self.__accuracyLimit


    ## Returns the maximal grid size
    # @return: maximal grid size
    def getGridSizeLimit(self):
        return self.__gridSizeLimit

    
    
    ## Checks if learning process have to be stopped
    # @param sampler: Learner object
    # @return: boolean value, true if learning has to stop, false otherwise
    def isTrainingComplete(self, sampler):
        ans = self.hasLimitReached(sampler) or not self.hasGridSizeChanged(sampler)
        self.__oldGridSize = sampler.getGrid().getSize()
        return ans

    
    def hasGridSizeChanged(self, sampler):
        return self.__oldGridSize != sampler.getGrid().getSize()
    
    
    def hasLimitReached(self, sampler):
        ans = (self.__adaptiveIterationLimit is None or
               self.__adaptiveIterationLimit < sampler.getCurrentIterationNumber()) \
               and (self.getGridSizeLimit() is None or
                    self.getGridSizeLimit() <= sampler.getGrid().getSize())
        return ans


    # # Setter for Maximal number of refinement iterations
    # @param limit: integer Maximal number of refinement iterations
    def setAdaptiveIterationLimit(self, limit):
        self.__adaptiveIterationLimit = limit


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
    
    
    # # Restores the SamplerStopPolicy object from the json object with attributes.
    #
    # @param cls python keyword (do not specify)
    # @param jsonObject A json object.
    # @return The restored SamplerStopPolicy object.
    @classmethod
    def fromJson(cls, jsonObject):
        policy = SamplerStopPolicy()
        if jsonObject.isContaining('_SamplerStopPolicy__adaptiveIterationLimit'):
            policy.setAdaptiveIterationLimit(jsonObject['_SamplerStopPolicy__adaptiveIterationLimit'])
        if jsonObject.isContaining('_SamplerStopPolicy__accuracyLimit'):
            policy.setAccuracyLimit(jsonObject['_SamplerStopPolicy__accuracyLimit'])
        if jsonObject.isContaining('_SamplerStopPolicy__gridSize'):
            policy.setGridSizeLimit(jsonObject['_SamplerStopPolicy__gridSize'])
        if jsonObject.isContaining('_SamplerStopPolicy__oldGridSize'):
            policy.__oldGridSize = jsonObject['_SamplerStopPolicy__oldGridSize']
        return policy

            
        

