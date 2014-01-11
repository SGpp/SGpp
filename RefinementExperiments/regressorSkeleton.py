#!/usr/bin/python

#set path
import sys
sys.path.append('/home/michael/workspace/SGpp')

# import modules
import bin.pysgpp
from bin.learner import LearnerBuilder, Regressor

import unittest
from operator import itemgetter
from numpy import zeros
from bin.learner.formatter.GridImageFormatter import GridImageFormatter
from bin.learner.formatter.GridFormatter import GridFormatter
from bin.pysgpp import DataMatrix, DataVector, RefinementDecorator
from bin.learner.solver.CGSolver import CGSolver
import pickle
from bin.pysgpp import DMSystemMatrix
import numpy as np
#from pylab import pcolor, colorbar, savefig, clf, get_cmap, close


from bin.pysgpp import HashRefinementBoundaries, RefinementDecorator, Grid, \
     createOperationMultipleEval, SurplusRefinementFunctor, SurplusCoarseningFunctor, \
     HashRefinement,SurplusVolumeRefinementFunctor, ANOVARefinement
from bin.learner.LearnerBuilder import LearnerBuilder
from bin.learner import Types,  LearnerEvents
from bin.data.DataContainer import DataContainer
import __main__


class RefinementTest:
    
    def setUp(self):
        self.example = "sum"
        self.refinement_functor = SurplusVolumeRefinementFunctor
        self.coarsening_functor = SurplusCoarseningFunctor
    
    def __init__(self):
        self.example = "sum"
        self.refinement_functor = SurplusVolumeRefinementFunctor
        self.coarsening_functor = SurplusCoarseningFunctor

    def testANOVAS(self, suffix):
#        from bin.controller.InfoToScreen import InfoToScreen
        builder = LearnerBuilder()
        builder = builder.buildRegressor()
        learner = builder.withTrainingDataFromCSVFile('refinement_strategy_%s.csv.gz'%suffix)\
        .withGrid().withLevel(3)\
        .withBorder(Types.BorderTypes.NONE)\
        .withSpecification().withIdentityOperator().withAdaptThreshold(0.001)\
        .withAdaptRate(1.0)\
        .withLambda(0.0001)\
        .withCGSolver().withImax(500)\
        .withStopPolicy().withAdaptiveItarationLimit(20)\
        .andGetResult()
        
        
        learner.specification.setBOperator(createOperationMultipleEval(learner.grid,
                    learner.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))
        
        while True: #repeat until policy says "stop"
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_STARTED)

            #learning step
            learner.alpha = learner.doLearningIteration(learner.dataContainer)
            learner.knowledge.update(learner.alpha)

            #compress grid
            if learner.iteration == 0:
                generator = learner.grid.createGridGenerator()
                functor = self.coarsening_functor(
                          learner.alpha,
                          generator.getNumberOfRemovablePoints(),
                          learner.specification.getAdaptThreshold())
                generator.coarsen(functor, learner.alpha)
            
            #self.plotGrid(learner, suffix)
            
            storage = learner.grid.getStorage()

            #self.plot_grid_historgram(suffix, learner, storage)
            
            
            
            #formatter = GridImageFormatter()
            #formatter.serializeToFile(learner.grid, "%s%d_projections_anova.png"%(suffix, learner.iteration))
            
            
            #calculate avg. error for training and test data and avg. for refine alpha
            learner.updateResults(learner.alpha, learner.dataContainer)
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            p_val = learner.trainAccuracy[-1] + learner.specification.getL()*np.sum(learner.alpha.array()**2)
            print "ANOVA %s iteration %d: %d grid points, %1.9f MSE, p* = %1.10f" % \
            (suffix, learner.iteration, storage.size(), learner.trainAccuracy[-1], p_val)           
            learner.iteration += 1
            if learner.iteration == 5: 
                pass
            if(learner.stopPolicy.isTrainingComplete(learner)): break
            
            #refine grid
            learner.notifyEventControllers(LearnerEvents.REFINING_GRID)
            learner.grid.getStorage().recalcLeafProperty()
            refinable_poits = learner.grid.createGridGenerator().getNumberOfRefinablePoints()
            pointsNum = learner.specification.getNumOfPointsToRefine(refinable_poits)
#            learner.grid.createGridGenerator().refine( SurplusRefinementFunctor(learner.errors, int(pointsNum), learner.specification.getAdaptThreshold()) )

            
            refiner = HashRefinement()
            functor = self.refinement_functor(learner.alpha, int(pointsNum), learner.specification.getAdaptThreshold())
            anova_refinement = ANOVARefinement(refiner)
            anova_refinement.free_refine(learner.grid.getStorage(), functor)
        #formatter = GridFormatter()
        #formatter.serializeToFile(learner.grid, "grid_anova_%s.txt"%suffix)

        del learner
        
    def testHashRefinement(self,suffix):
#        from bin.controller.InfoToScreen import InfoToScreen
        builder = LearnerBuilder()
        builder = builder.buildRegressor()
        learner = builder.withTrainingDataFromCSVFile('refinement_strategy_%s.csv.gz'%suffix)\
        .withGrid().withLevel(3)\
        .withBorder(Types.BorderTypes.NONE)\
        .withSpecification().withIdentityOperator().withAdaptThreshold(0.001)\
        .withAdaptRate(0.8)\
        .withLambda(0.00001)\
        .withCGSolver().withImax(10000000)\
        .withStopPolicy().withAdaptiveItarationLimit(20)\
        .andGetResult()
        
        
        learner.specification.setBOperator(createOperationMultipleEval(learner.grid,
                    learner.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))
        
        while True: #repeat until policy says "stop"
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_STARTED)

            #learning step
            learner.alpha = learner.doLearningIteration(learner.dataContainer)
            learner.knowledge.update(learner.alpha)

            # #compress grid
            # if learner.iteration == 0:
            #     generator = learner.grid.createGridGenerator()
            #     functor = self.coarsening_functor(
            #               learner.alpha,
            #               generator.getNumberOfRemovablePoints(),
            #               learner.specification.getAdaptThreshold())
            #     generator.coarsen(functor, learner.alpha)
            
            storage = learner.grid.getStorage()
            
            #calculate avg. error for training and test data and avg. for refine alpha
            learner.updateResults(learner.alpha, learner.dataContainer)
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            p_val = learner.trainAccuracy[-1] + learner.specification.getL()*np.sum(learner.alpha.array()**2)
            print "HashRefinement %s iteration %d: %d grid points, %1.9f MSE, p* = %1.10f" % \
            (suffix, learner.iteration, storage.size(), learner.trainAccuracy[-1], p_val)           
            learner.iteration += 1
            if learner.iteration == 5: 
                pass
            if(learner.stopPolicy.isTrainingComplete(learner)): break
            
            #refine grid
            learner.notifyEventControllers(LearnerEvents.REFINING_GRID)
            learner.grid.getStorage().recalcLeafProperty()
            refinable_poits = learner.grid.createGridGenerator().getNumberOfRefinablePoints()
            pointsNum = learner.specification.getNumOfPointsToRefine(refinable_poits)
            
            refiner = HashRefinement()
            self.refinement_functor = SurplusRefinementFunctor
            functor = self.refinement_functor(learner.alpha, int(pointsNum), learner.specification.getAdaptThreshold())
            #anova_refinement = ANOVARefinement(refiner)
            #anova_refinement.free_refine(learner.grid.getStorage(), functor)
            refiner.free_refine(learner.grid.getStorage(), functor)

        del learner

def main():
    test = RefinementTest()
    test.testANOVAS('sum')
    print
    test.testHashRefinement('sum')
    
    
if __name__ == "__main__":
    main()

  



