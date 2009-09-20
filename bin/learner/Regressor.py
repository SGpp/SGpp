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
# @brief Regression
# @version $HEAD$

from Learner import Learner, LearnerEvents
from bin.pysgpp import DataVector, SurplusRefinementFunctor

class Regressor(Learner):
    
    errors = None # Errors per basis function
    
    ## Evaluate regression MSE
    #
    # @param data: DataContainer dataset
    # @param alpha: DataVector alpha-vector
    # @return: mean square error
    def evalError(self, data, alpha):
        size = data.getPoints().getSize()
        if size == 0: return 0
        
        error = DataVector(size)
        self.specification.getBOperator().multTranspose(alpha, data.getPoints(), error)
        error.sub(data.getValues()) # error vector
        error.sqr() # entries squared
        mse = error.sum() / size # MSE
        
        # calculate error per basis function
        if self.errors == None:
            self.errors = DataVector(size)
        self.specification.getBOperator().mult(error, data.getPoints(), self.errors)
        
        return mse
    
    
    ##Update different statistics about training progress
    # @param alpha: DataVector alpha-vector
    # @param trainSubset: DataContainer with training data
    # @param testSubset: DataContainer with validation data, default value: None
    def updateResults(self, alpha, trainSubset, testSubset = None):
        self.knowledge.update(alpha)
        # @todo (khakhutv) Add L2-norm of error as well as min/max errors
        #eval Error for training data and append it to other in this iteration
        self.trainAccuracy.append(self.evalError(trainSubset, alpha))
        
        i = float(len(self.trainAccuracy))
        
        #eval error for test data and append it to other in this iteration
        if testSubset != None:  
            self.testAccuracy.append(self.evalError(testSubset, alpha))
            self.testingOverall.append(sum(self.testAccuracy)/i)
            
        self.trainingOverall.append(sum(self.trainAccuracy)/i)

        # @todo (khakhutv) grid.getSize() change Grid interface
        self.numberPoints.append(self.grid.getStorage().size())
    
    
    ##Refines grid with the number of points as specified in corresponding TrainingSpecification object
    def refineGrid(self):
        self.notifyEventControllers(LearnerEvents.REFINING_GRID)
        
        pointsNum = self.specification.getNumOfPointsToRefine( self.grid.createGridGenerator().getNumberOfRefinablePoints() )
        
        # @todo (khakhutv) develop a way to simplify interfaces and use different functors
        self.grid.createGridGenerator().refine( SurplusRefinementFunctor(self.errors, pointsNum, self.specification.getAdaptThreshold()) )
        