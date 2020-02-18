# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

###############################################################################
# $Copyright$ #
###############################################################################

import unittest
from operator import itemgetter
from numpy import zeros
from bin.learner.formatter.GridImageFormatter import GridImageFormatter
from bin.learner.formatter.GridFormatter import GridFormatter
from pysgpp import DataMatrix, DataVector
from bin.learner.solver.CGSolver import CGSolver
import pickle
from pysgpp import DMSystemMatrix
import numpy as np
#from pylab import pcolor, colorbar, savefig, clf, get_cmap, close


from pysgpp import HashRefinementBoundaries, ANOVARefinement, Grid, \
     createOperationMultipleEval, SurplusRefinementFunctor, SurplusCoarseningFunctor, \
     HashRefinement,SurplusVolumeRefinementFunctor, ANOVACoarseningFunctor
from bin.learner.LearnerBuilder import LearnerBuilder
from bin.learner import Types,  LearnerEvents
from bin.data.DataContainer import DataContainer


class TestRefinementCoerseningANOVAStrategy(unittest.TestCase):
    
    def setUp(self):
        self.example = "sum"
        self.refinement_functor = SurplusVolumeRefinementFunctor
        self.coarsening_functor = ANOVACoarseningFunctor #SurplusCoarseningFunctor
        
#    def tearDown(self):
#        import gc
#        gc.collect()
#        plt.close('all')


    def plot_grid_historgram(self, suffix, learner, storage, exp_type='anova'):
        import matplotlib.pyplot as plt

        m = self.process_grid_statistics(storage)
        figure = plt.figure()
        plt.pcolor(m, cmap=plt.get_cmap('RdYlGn'), figure=figure)
        plt.colorbar()
        plt.savefig("%s%d_%s.png" % (suffix, learner.iteration, exp_type), figure=figure)
        plt.close(figure)

        
    def plotGrid(self, learner, suffix):
        from mpl_toolkits.mplot3d.axes3d import Axes3D
        import matplotlib.pyplot as plt
        xs = np.linspace(0, 1, 30)
        ys = np.linspace(0, 1, 30)
        X, Y = np.meshgrid(xs, ys)
        Z = zeros(np.shape(X))
        input = DataMatrix(np.shape(Z)[0]*np.shape(Z)[1], 2)
        r = 0
        for i in range(np.shape(Z)[0]):
            for j in range(np.shape(Z)[1]):
                input.set(r, 0, X[i,j])
                input.set(r, 1, Y[i,j])
                r += 1
        result = learner.applyData(input)
        r = 0
        for i in range(np.shape(Z)[0]):
            for j in range(np.shape(Z)[1]):
                Z[i,j] = result[r]
                r += 1
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot_wireframe(X,Y,Z) 
        #plt.draw()
        plt.savefig("grid3d_%s_%i.png" % (suffix, learner.iteration))
        fig.clf()
        plt.close(plt.gcf())


    def testANOVASum(self):
        self.__testANOVAS('sum')

    def _testANOVAX1(self):
        self.__testANOVAS('x1')

    def _testANOVASinCos(self):
        self.__testANOVAS('sincos')


    def _testSpaceSum(self):
        self.__testSpaceS('sum')

    def _testSpaceX1(self):
        self.__testSpaceS('x1')

    def _testSpaceSinCos(self):
        self.__testSpaceS('sincos')


        
    def process_grid_statistics(self, storage):
        grid_statistics = {}
        for i in range(storage.getSize()):
            point = storage.getPoint(i)
            key = (point.getLevel(0),point.getLevel(1))
            if key in grid_statistics:
                grid_statistics[key] += 1
            else:
                grid_statistics[key] = 1
        keys = sorted(grid_statistics.keys())
        
        l1max = max(keys, key=itemgetter(0))[0]
        l2max = max(keys, key=itemgetter(1))[1]
        m = zeros([l1max, l2max])
            
        for k in keys:
            m[k[0]-1, k[1]-1] = grid_statistics[k]
        return m
     


    def __testANOVAS(self, suffix):
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
        .withStopPolicy().withAdaptiveItarationLimit(5)\
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
                generator = learner.grid.getGenerator()
                functor = self.coarsening_functor(
                          learner.alpha,
                          generator.getNumberOfRemovablePoints(),
                          0.99, learner.grid.getStorage())
#                functor = self.coarsening_functor(
#                          learner.alpha,
#                          generator.getNumberOfRemovablePoints(),
#                          learner.specification.getAdaptThreshold())
                generator.coarsen(functor)
            #print "coersening finished"
            self.plotGrid(learner, suffix)
            
            storage = learner.grid.getStorage()

            self.plot_grid_historgram(suffix, learner, storage)
            
            
            
            formatter = GridImageFormatter()
            formatter.serializeToFile(learner.grid, "%s%d_projections_anova.png"%(suffix, learner.iteration))
            
            
            #calculate avg. error for training and test data and avg. for refine alpha
            learner.updateResults(learner.alpha, learner.dataContainer)
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            p_val = learner.trainAccuracy[-1] + learner.specification.getL()*np.sum(learner.alpha.array()**2)
            print("ANOVA %s iteration %d: %d grid points, %1.9f MSE, p* = %1.10f" % \
            (suffix, learner.iteration, storage.getSize(), learner.trainAccuracy[-1], p_val))           
            learner.iteration += 1
            if learner.iteration == 5: 
                pass
            if(learner.stopPolicy.isTrainingComplete(learner)): break
            
            #refine grid
            learner.notifyEventControllers(LearnerEvents.REFINING_GRID)
            learner.grid.getStorage().recalcLeafProperty()
            refinable_poits = learner.grid.getGenerator().getNumberOfRefinablePoints()
            pointsNum = learner.specification.getNumOfPointsToRefine(refinable_poits)
#            learner.grid.getGenerator().refine( SurplusRefinementFunctor(learner.errors, int(pointsNum), learner.specification.getAdaptThreshold()) )

            
            refiner = HashRefinement()
            functor = self.refinement_functor(learner.alpha, int(pointsNum), learner.specification.getAdaptThreshold())
            anova_refinement = ANOVARefinement(refiner)
            anova_refinement.free_refine(learner.grid.getStorage(), functor)
        #formatter = GridFormatter()
        #formatter.serializeToFile(learner.grid, "grid_anova_%s.txt"%suffix)

        del learner


    def __testSpaceS(self, suffix):
        
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
        .withStopPolicy().withAdaptiveItarationLimit(5)\
        .andGetResult()
        
        
        learner.specification.setBOperator(createOperationMultipleEval(learner.grid,
                    learner.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))
        
        while True: #repeat until policy says "stop"
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_STARTED)

            #learning step
            learner.alpha = learner.doLearningIteration(learner.dataContainer)
            learner.knowledge.update(learner.alpha)
            
            #self.plotGrid(learner, suffix)
            
            storage = learner.grid.getStorage()

            #self.plot_grid_historgram(suffix, learner, storage, 'space')
            
            
            
#            formatter = GridImageFormatter()
#            formatter.serializeToFile(learner.grid, "%s%d_projections_space.png"%(suffix, learner.iteration))
            
            
            #calculate avg. error for training and test data and avg. for refine alpha
            learner.updateResults(learner.alpha, learner.dataContainer)
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            p_val = learner.trainAccuracy[-1] + learner.specification.getL()*np.sum(learner.alpha.array()**2)
            print("Space %s iteration %d: %d grid points, %1.9f MSE, p* = %1.10f" % \
            (suffix, learner.iteration, storage.getSize(), learner.trainAccuracy[-1], p_val))           
            learner.iteration += 1
            if(learner.stopPolicy.isTrainingComplete(learner)): break
            
            #refine grid
            learner.notifyEventControllers(LearnerEvents.REFINING_GRID)
        
            pointsNum = learner.specification.getNumOfPointsToRefine( learner.grid.getGenerator().getNumberOfRefinablePoints() )
            learner.grid.getGenerator().refine( self.refinement_functor(learner.alpha, int(pointsNum), learner.specification.getAdaptThreshold()) )
        #formatter = GridFormatter()
        #formatter.serializeToFile(learner.grid, "grid_anova_%s.txt"%suffix)

        del learner
    
    
    
           
        
        
        
        
        
        
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()
