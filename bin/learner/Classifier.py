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

## @package Learner
# @ingroup bin.learner
# @brief Classifications
# @version $CURR$

from Learner import Learner, LearnerEvents


class Classifier(Learner):


    ## Evaluate classification accuracy as percent of correct classified data points
    #
    # @param data: DataContainer dataset
    # @param alpha: DataVector alpha-vector
    # @return: percent of correct classified data, 0 if data set is empty
    def evalError(self, data, alpha):
        evalOp = self.grid.createOperationEval()
        size = data.getPoints().getSize()
        if size != 0: return evalOp.test(alpha, data.getPoints(), data.getValues()) / size
        else: return 0
    
    
    ##Update different statistics about training progress
    #
    # @param alpha: DataVector alpha-vector
    # @param trainSubset: DataContainer with training data
    # @param testSubset: DataContainer with validation data, default value: None
    def updateResults(self, alpha, trainSubset, testSubset = None):
        #eval Error for training data and append it to other in this iteration
        self.trainAccuracy.append(self.evalError(trainSubset, alpha))
        
        i = float(len(self.trainAccuracy))
        
        #eval error for test data and append it to other in this iteration
        if testSubset != None:  
            self.testAccuracy.append(self.evalError(testSubset, alpha))
            self.testingOverall.append(sum(self.testAccuracy)/i)
            
        self.trainingOverall.append(sum(self.trainAccuracy)/i)

        #@todo: grid.getSize() change Grid interface
        self.numberPoints.append(self.grid.getStorage().size())
    
    
    ##Refines grid with the number of points as specified in corresponding TrainingSpecification object
    def refineGrid(self):
        self.notifyEventControllers(LearnerEvents.REFINING_GRID)
        self.grid.refine(self.alpha, 
                             self.specification.getNumOfPointsToRefine(
                                          self.grid.createGridGenerator().getNumberOfRefinablePoints()
                                                                      )
                             )


