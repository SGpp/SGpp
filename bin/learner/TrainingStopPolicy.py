##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2009 Valeriy Khakhutskyy (khakhutv@in.tum.de)               #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################



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
        if (self.__adaptiveIterationLimit == None 
        or self.__adaptiveIterationLimit < learner.getCurrentIterationNumber()) \
        and (self.getGridSizeLimit() == None 
        or self.getGridSizeLimit() <= learner.grid.getSize()) \
        and (self.getMSELimit() == None
        or self.getMSELimit() >= learner.trainingOverall[-1]):
            return True
        if learner.grid.getSize() == self.__oldGridSize: return True
        self.__oldGridSize = learner.grid.getSize()
        return False

    
    ## Setter for Maximal number of refinement iterations
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

            
        


