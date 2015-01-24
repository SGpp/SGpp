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


from Learner import Learner, LearnerEvents
from pysgpp import SurplusRefinementFunctor, createOperationTest


## The class implements the abstract methods from Learner and allows to accomplish
# basic classification tasks.
class Classifier(Learner):


    ## Evaluate classification accuracy as percent of correct classified data points
    #
    # @param data: DataContainer dataset
    # @param alpha: DataVector alpha-vector
    # @return: percent of correct classified data, 0 if data set is empty
    def evalError(self, data, alpha):
        testOp = createOperationTest(self.grid)
        size = data.getPoints().getNrows()
        if size != 0: return testOp.test(alpha, data.getPoints(), data.getValues()) / size
        else: return 0
    
    
    ##Update different statistics about training progress
    #
    # @param alpha: DataVector alpha-vector
    # @param trainSubset: DataContainer with training data
    # @param testSubset: DataContainer with validation data, default value: None
    def updateResults(self, alpha, trainSubset, testSubset = None):
        self.knowledge.update(alpha)
        #eval Error for training data and append it to other in this iteration
        self.trainAccuracy.append(self.evalError(trainSubset, alpha))
        
        i = float(len(self.trainAccuracy))
        
        #eval error for test data and append it to other in this iteration
        if testSubset != None:  
            self.testAccuracy.append(self.evalError(testSubset, alpha))
            self.testingOverall.append(sum(self.testAccuracy)/i)
            
        self.trainingOverall.append(sum(self.trainAccuracy)/i)

        self.numberPoints.append(self.grid.getSize())
    
    
    ##Refines grid with the number of points as specified in corresponding TrainingSpecification object
    def refineGrid(self):
        self.notifyEventControllers(LearnerEvents.REFINING_GRID)
        pointsNum = self.specification.getNumOfPointsToRefine( self.grid.createGridGenerator().getNumberOfRefinablePoints() )
        
        # @todo (khakhutv) (low) develop a way to simplify interfaces and use different functors
        self.grid.createGridGenerator().refine( SurplusRefinementFunctor(self.alpha, pointsNum, self.specification.getAdaptThreshold()) )
        
    ## Creates Classifier object without any parameter set
    # @param jsonObject: A json object.
    # @return: Classifier object without any parameter set
#    @classmethod
#    def fromJson(cls, jsonObject):
#        return Classifier()
                            


