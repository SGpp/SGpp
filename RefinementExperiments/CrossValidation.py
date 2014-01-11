#!/usr/bin/python

#set path
import sys
import os
sys.path.append('/home/michael/workspace/SG++')

import bin.pysgpp
from bin.learner import LearnerBuilder, Regressor

import unittest
from operator import itemgetter
from numpy import zeros
from bin.learner.formatter.GridImageFormatter import GridImageFormatter
from bin.learner.formatter.GridFormatter import GridFormatter
from bin.pysgpp import DataMatrix, DataVector, RefinementDecorator
from bin.controller.InfoToScreen import InfoToScreen
from bin.controller.InfoToScreenRegressor import InfoToScreenRegressor
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

class CrossValidationTest:
    """docstring for CrossValidationTest"""
    
    def __init__(self, kwargs):
        self.kwargs = kwargs

    #
    #grid point based
    #

    #hash refinement
    def testHashRefinement(self, kwargs):
        regressor = LearnerBuilder()
        regressor = regressor.buildRegressor()\
                .withTrainingDataFromARFFFile(self.kwargs['training'])\
                .withTestingDataFromARFFFile(self.kwargs['testing'])\
                .withGrid()\
                    .withLevel(self.kwargs['gridLevel'])\
                    .withBorder(Types.BorderTypes.NONE)\
                .withSpecification()\
                    .withIdentityOperator()\
                    .withAdaptPoints(kwargs['refinementsNum'])\
                    .withLambda(kwargs['lambda'])\
                .withCGSolver()\
                    .withImax(500)\
                .withStopPolicy()\
                    .withAdaptiveItarationLimit(self.kwargs['maxRefinements'])\
                .withProgressPresenter(InfoToScreenRegressor())\
                .andGetResult()

        regressor.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STARTED)
        regressor.specification.setBOperator(createOperationMultipleEval(regressor.grid,
                                        regressor.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))

        #get corresponding datasets
        trainSubset = regressor.dataContainer.getTrainDataset()
        testSubset = regressor.dataContainer.getTestDataset()

        while True:
            #start learning
            regressor.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STEP_STARTED)

            regressor.alpha = regressor.doLearningIteration(trainSubset)

            #calculate avg. error for training and test data and avg. for refine alpha
            regressor.updateResults(regressor.alpha, trainSubset, testSubset)

            regressor.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE)
            regressor.iteration += 1

            if(regressor.stopPolicy.isTrainingComplete(regressor)): break

            #refine grid
            regressor.notifyEventControllers(LearnerEvents.REFINING_GRID)
            regressor.grid.getStorage().recalcLeafProperty()

            pointsNum = regressor.specification.getAdaptPoints()
            
            refiner = HashRefinement()
            functor = SurplusRefinementFunctor(regressor.alpha, int(pointsNum), 0)
            refiner.free_refine(regressor.grid.getStorage(), functor)

        del regressor
        
        
    #predictive refinement
    
    #anova refinement

    #predictive anova refinement

    #predictive strack anova refinement

    #
    #subspace based refinement
    #

    #subspace refinement

    #subspace gsg refinement

    #predictive subspace gsg refinement


def main():

    #create argumentlist:
    commonKWArgs = {'training': 'friedman2_10000_train.arff.gz' ,
    'testing': 'friedman2_10000_test.arff.gz',
    'gridLevel': 3,
    'maxIterCG': 1000,
    'maxRefinements': 4
    }
    standardPointKWArgs = {'refinementsNum': 50,
                            'lambda': 0.0001,
                            }

    tests = CrossValidationTest(commonKWArgs)
    tests.testHashRefinement(standardPointKWArgs)



if __name__ == "__main__":
    main()
