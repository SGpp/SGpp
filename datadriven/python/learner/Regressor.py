# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


from Learner import Learner, LearnerEvents
from pysgpp import DataVector, SurplusRefinementFunctor

from  math import sqrt


## Subclass of Learner, responsible for regression.
# The methods specific for regression are implemented here.
class Regressor(Learner):
    
    ## Errors per basis function
    errors = None 
    
    ## Error vector
    error = None
    
    
    ##constructor
    def __init__(self):
        super(Regressor,self).__init__()
       
        
    ##calculate L2-norm of error
    # @return: last L2-norm of error
    def getL2NormError(self):
        return sqrt(self.error.sum())
    
    
    ## calculate max error
    # @return: max error
    def getMaxError(self):
        return sqrt(self.error.max())
    
    
    ## calculate min error
    # @return: min error
    def getMinError(self):
        return sqrt(self.error.min())
    
    
    ## Evaluate regression MSE
    #
    # @param data: DataContainer dataset
    # @param alpha: DataVector alpha-vector
    # @return: mean square error
    def evalError(self, data, alpha):
        size = data.getPoints().getNrows()
        if size == 0: return 0
        
        self.error = DataVector(size)
        self.specification.getBOperator(data.getName()).mult(alpha, self.error)
        self.error.sub(data.getValues()) # error vector
        self.error.sqr() # entries squared
        errorsum = self.error.sum()
        mse = errorsum / size # MSE
        
        # calculate error per basis function
        self.errors = DataVector(len(alpha))
        self.specification.getBOperator(data.getName()).multTranspose(self.error, self.errors)
        self.errors.componentwise_mult(alpha)
        
        return mse
    
    
    ##Update different statistics about training progress
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
        
        pointsNum = self.specification.getNumOfPointsToRefine( self.grid.getGenerator().getNumberOfRefinablePoints() )
        self.grid.getGenerator().refine( SurplusRefinementFunctor(self.errors, pointsNum, self.specification.getAdaptThreshold()) )

