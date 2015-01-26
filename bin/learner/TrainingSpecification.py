# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

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

## Collection of parameters, which specify the learning process.
# An object of the class is aggregated by the Learner object.
class TrainingSpecification(object):

    __adaptPoints = 0       #Number of points to refine in one refinement iteration
    __l = None              #Regularization parameter
    __adaptRate = 0         #Rate of points to refine in one refinement iteration, between 0 and 1
    __adaptThreshold = 0.0   #threshold, only the points with greater to equal absolute values of the refinement criterion (e.g. alpha or error) will be refined
    __cOperator = None      #C operator
    __bOperator = None      #B operator
    __cOperatorType = None  #Type of the c operator as a string
    __vecType = None
    
    ## Returns the type of the C operator
    # @return: the type of the C operator as a string
    def getCOperatorType(self):
        return self.__cOperatorType

    ## Sets the type of the C operator
    # @param value string type of the C operator
    def setCOperatorType(self, value):
        self.__cOperatorType = value

    
    ## Getter for refinement threshold
    # only the points with greater to equal absolute values of 
    # the refinement criterion (e.g. alpha or error) will be refined
    #
    # @return float threshold
    def getAdaptThreshold(self):
        return self.__adaptThreshold

    ## Setter for refinement threshold
    # only the points with greater to equal absolute values of 
    # the refinement criterion (e.g. alpha or error) will be refined
    #
    # @param value: float threshold
    def setAdaptThreshold(self, value):
        self.__adaptThreshold = value
    
    ## Setter for Number of points to refine
    #
    # @param value: integer Number of points to refine
    def setAdaptPoints(self, value):
        self.__adaptPoints = value


    ## Setter for Regularization parameter
    #
    # @param value: double Regularization parameter
    def setL(self, value):
        self.__l = value


    ## Setter for Rate of points to refine
    #
    # @param value: double in [0,1] Rate of points to refine
    def setAdaptRate(self, value):
        self.__adaptRate = value


    ## Setter for C operator
    #
    # @param value: OperationMatrix
    def setCOperator(self, value):
        self.__cOperator = value


    ## Setter for B operator
    #
    # @param value: OperationB
    # @param name: operator identifier
    def setBOperator(self, value, name="train"):
        if self.__bOperator == None: self.__bOperator = {}
        self.__bOperator[name] = value


    ## Getter for Number of points to refine
    #
    # @return: integer Number of points to refine
    def getAdaptPoints(self):
        return self.__adaptPoints


    ## Getter for Regularization parameter
    #
    # @return: double Regularization parameter
    def getL(self):
        return self.__l


    ## Getter for Rate of points to refine
    #
    # @return: double in [0,1] Rate of points to refine
    def getAdaptRate(self):
        return self.__adaptRate


    ## Getter for C operator
    #
    # @return: OperationMatrix
    def getCOperator(self):
        return self.__cOperator


    ## Getter for B operator
    #
    # @param name: operator identifier
    # @return: OperationB
    def getBOperator(self, name="train"):
        if self.__bOperator != None:
            return self.__bOperator[name]
        else: return None
        
        
    def getVectorizationType(self):
        return self.__vecType
    
    
    def setVectorizationType(self, vecType):
        self.__vecType = vecType
    
    
    ## Calculates the number of points which should be refined
    #
    # @param refinablePoints: integer number of points which can be refined
    # @return: integer number of point which should be refined
    def getNumOfPointsToRefine(self, refinablePoints):
        ratePoints = self.__adaptRate * refinablePoints
        if self.__adaptPoints == 0:
            return ratePoints
        elif self.__adaptRate == 0:
            return self.__adaptPoints
        else: 
            return min(ratePoints, self.__adaptPoints)
        
        
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object.   
    def toString(self):
        serializationString = "'module' : '" + self.__module__ + "',\n"
        for attrName in dir(self):
            attrValue = self.__getattribute__(attrName)
            
            #lists, integers, floats, dictionaries can serialized with str()
            if type(attrValue) in [types.ListType, types.IntType, 
                             types.FloatType] and attrName.find("__") != 0: 
                serializationString += "'" + attrName + "'" + " : " + str(attrValue) + ",\n"
                
            # serialize strings with quotes    
            elif type(attrValue) == types.StringType and attrName.find("__") != 0:
                serializationString += "'" + attrName + "'" + " : '" + attrValue + "',\n"

        serializationString = "{" + serializationString.rstrip(",\n") + "}"
        return serializationString
    
    
    ## Restores the TrainingSpecification object from the json object with attributes.
    #
    # @param cls python keyword (do not specify)
    # @param jsonObject A json object.
    # @return The restored TrainingSpecification object.
    @classmethod
    def fromJson(cls, jsonObject):
        specification = TrainingSpecification()
        if jsonObject.has_key('_TrainingSpecification__adaptPoints'):
            specification.__adaptPoints = jsonObject['_TrainingSpecification__adaptPoints']
        if jsonObject.has_key('_TrainingSpecification__l'):
            specification.__l = jsonObject['_TrainingSpecification__l']
        if jsonObject.has_key('_TrainingSpecification__adaptRate'):
            specification.__adaptRate = jsonObject['_TrainingSpecification__adaptRate']
        if jsonObject.has_key('_TrainingSpecification__adaptThreshold'):
            specification.__adaptThreshold = jsonObject['_TrainingSpecification__adaptThreshold']
        if jsonObject.has_key('_TrainingSpecification__vecType'):
            specification.__vecType = jsonObject['_TrainingSpecification__vecType']
        specification.__cOperator = None
        specification.__bOperator = None
        return specification




