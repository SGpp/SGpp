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

## @package TrainingStopPolicy
# @ingroup bin.learner
# @brief Rules to stop learning iterations to prevent overfitting
# @version $CURR$

class TrainingStopPolicy(object):
    __adaptiveIterationLimit = None     #Maximal number of refinement iterations
    __epochsLimit = None                #Maximal number of iterations, during which accuracy can decreases
    __MSELimit = None                   #MSE on validation data, that have to be achieved
    __accuracyLimit = None              #accuracy on validation data, that have to be achieved
    __gridSize = None                   #Maximal grid size
    
    
    ##Contructor
    def __init__(self):
        self.__adaptiveIterationLimit = None
        self.__epochsLimit = None
        self.__MSELimit = None
        self.__accuracyLimit = None
        self.__gridSize = None
    
    
    ## Checks if learning process have to be stopped
    # @todo: Make it more advanced
    # @param learner: Learner object 
    # @return: boolean value, true if learning has to stop, false otherwise
    def isTrainingComplete(self, learner):
        if self.__adaptiveIterationLimit != None and self.__adaptiveIterationLimit < learner.getCurrentIterationNumber():
            return True
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
        self.__gridSize = limit
        


