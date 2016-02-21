# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

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
        pointsNum = self.specification.getNumOfPointsToRefine( self.grid.getGenerator().getNumberOfRefinablePoints() )
        
        self.grid.getGenerator().refine( SurplusRefinementFunctor(self.alpha, pointsNum, self.specification.getAdaptThreshold()) )
        

