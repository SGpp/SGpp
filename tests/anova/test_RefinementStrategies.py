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
from bin.pysgpp import DMSystemMatrix
import numpy as np
#from pylab import pcolor, colorbar, savefig, clf, get_cmap, close

from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
        
class TestRefinementANOVAStrategy(unittest.TestCase):
    
    def setUp(self):
        self.example = "sum"
        
#    def tearDown(self):
#        import gc
#        gc.collect()
#        plt.close('all')
        
    def plotGrid(self, learner, suffix):
#        plt.ioff()
        xs = np.linspace(0, 1, 30)
        ys = np.linspace(0, 1, 30)
        X, Y = np.meshgrid(xs, ys)
        Z = zeros(np.shape(X))
        input = DataMatrix(np.shape(Z)[0]*np.shape(Z)[1], 2)
        r = 0
        for i in xrange(np.shape(Z)[0]):
            for j in xrange(np.shape(Z)[1]):
                input.set(r, 0, X[i,j])
                input.set(r, 1, Y[i,j])
                r += 1
        result = learner.applyData(input)
        r = 0
        for i in xrange(np.shape(Z)[0]):
            for j in xrange(np.shape(Z)[1]):
                Z[i,j] = result[r]
                r += 1
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot_wireframe(X,Y,Z) 
        #plt.draw()
        plt.savefig("grid3d_%s_%i.png" % (suffix, learner.iteration))
        fig.clf()
        plt.close(plt.gcf())


    def _testANOVASum(self):
        self.__testANOVAS('sum')

    def _testANOVAX1(self):
        self.__testANOVAS('x1')

    def _testANOVASinCos(self):
        self.__testANOVAS('sincos')


    def _testSpaceSum(self):
        self.__testSpaceS('sum')

    def _testSpaceAX1(self):
        self.__testSpaceS('x1')

    def testSpaceSinCos(self):
        self.__testSpaceS('sincos')


        
    def process_grid_statistics(self, storage):
        grid_statistics = {}
        for i in xrange(storage.size()):
            point = storage.get(i)
            key = (point.getLevel(0),point.getLevel(1))
            if grid_statistics.has_key(key):
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
     

    def plot_grid_historgram(self, suffix, learner, storage, exp_type='anova'):
        m = self.process_grid_statistics(storage)
        figure = plt.figure()
        plt.pcolor(m, cmap=plt.get_cmap('RdYlGn'), figure=figure)
        plt.colorbar()
        plt.savefig("%s%d_%s.png" % (suffix, learner.iteration, exp_type), figure=figure)
        plt.close(figure)


    def __testANOVAS(self, suffix):
        from bin.pysgpp import HashRefinementBoundaries, RefinementANOVAStrategy, Grid, createOperationMultipleEval,\
        SurplusRefinementFunctor, HashRefinement, SurplusCoarseningFunctor
        from bin.learner.LearnerBuilder import LearnerBuilder
        from bin.learner import Types,  LearnerEvents
        from bin.data.DataContainer import DataContainer
#        from bin.controller.InfoToScreen import InfoToScreen
        builder = LearnerBuilder()
        builder = builder.buildRegressor()
        learner = builder.withTrainingDataFromCSVFile('refinement_strategy_%s.csv'%suffix)\
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
                generator = learner.grid.createGridGenerator()
                functor = SurplusCoarseningFunctor(
                          learner.alpha,
                          generator.getNumberOfRemovablePoints(),
                          learner.specification.getAdaptThreshold())
                generator.coarsen(functor, learner.alpha)
            
            self.plotGrid(learner, suffix)
            
            storage = learner.grid.getStorage()

            self.plot_grid_historgram(suffix, learner, storage)
            
            
            
            formatter = GridImageFormatter()
            formatter.serializeToFile(learner.grid, "%s%d_projections_anova.png"%(suffix, learner.iteration))
            
            
            #calculate avg. error for training and test data and avg. for refine alpha
            learner.updateResults(learner.alpha, learner.dataContainer)
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            p_val = learner.trainAccuracy[-1] + learner.specification.getL()*np.sum(learner.alpha.array()**2)
            print "ANOVA %s iteration %d: %d grid points, %1.9f MSE, p* = %1.10f" % \
            (suffix, learner.iteration, storage.size(), learner.trainAccuracy[-1], p_val)           
            learner.iteration += 1
            if(learner.stopPolicy.isTrainingComplete(learner)): break
            
            #refine grid
            learner.notifyEventControllers(LearnerEvents.REFINING_GRID)
        
            pointsNum = learner.specification.getNumOfPointsToRefine( learner.grid.createGridGenerator().getNumberOfRefinablePoints() )
#            learner.grid.createGridGenerator().refine( SurplusRefinementFunctor(learner.errors, int(pointsNum), learner.specification.getAdaptThreshold()) )

            refiner = HashRefinement()
            functor = SurplusRefinementFunctor(learner.alpha, int(pointsNum), learner.specification.getAdaptThreshold())
            refinement_strategy = RefinementANOVAStrategy(functor)
            refiner.strategy_refine(learner.grid.getStorage(), refinement_strategy)
        #formatter = GridFormatter()
        #formatter.serializeToFile(learner.grid, "grid_anova_%s.txt"%suffix)

        del learner


    def __testSpaceS(self, suffix):
        from bin.pysgpp import HashRefinementBoundaries, RefinementANOVAStrategy, Grid, createOperationMultipleEval,\
        SurplusRefinementFunctor, HashRefinement, SurplusCoarseningFunctor
        from bin.learner.LearnerBuilder import LearnerBuilder
        from bin.learner import Types,  LearnerEvents
        from bin.data.DataContainer import DataContainer
#        from bin.controller.InfoToScreen import InfoToScreen
        builder = LearnerBuilder()
        builder = builder.buildRegressor()
        learner = builder.withTrainingDataFromCSVFile('refinement_strategy_%s.csv'%suffix)\
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

            self.plot_grid_historgram(suffix, learner, storage, 'space')
            
            
            
            formatter = GridImageFormatter()
            formatter.serializeToFile(learner.grid, "%s%d_projections_space.png"%(suffix, learner.iteration))
            
            
            #calculate avg. error for training and test data and avg. for refine alpha
            learner.updateResults(learner.alpha, learner.dataContainer)
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            p_val = learner.trainAccuracy[-1] + learner.specification.getL()*np.sum(learner.alpha.array()**2)
            print "Space %s iteration %d: %d grid points, %1.9f MSE, p* = %1.10f" % \
            (suffix, learner.iteration, storage.size(), learner.trainAccuracy[-1], p_val)           
            learner.iteration += 1
            if(learner.stopPolicy.isTrainingComplete(learner)): break
            
            #refine grid
            learner.notifyEventControllers(LearnerEvents.REFINING_GRID)
        
            pointsNum = learner.specification.getNumOfPointsToRefine( learner.grid.createGridGenerator().getNumberOfRefinablePoints() )
            learner.grid.createGridGenerator().refine( SurplusRefinementFunctor(learner.alpha, int(pointsNum), learner.specification.getAdaptThreshold()) )
        #formatter = GridFormatter()
        #formatter.serializeToFile(learner.grid, "grid_anova_%s.txt"%suffix)

        del learner
    
    
    def __testRefineANOVA(self):
        return
        from bin.pysgpp import HashRefinementBoundaries, RefinementANOVAStrategy, Grid, createOperationMultipleEval,\
        SurplusRefinementFunctor, HashRefinement, SurplusCoarseningFunctor
        from bin.learner.LearnerBuilder import LearnerBuilder
        from bin.learner import Types,  LearnerEvents
        from bin.data.DataContainer import DataContainer
#        from bin.controller.InfoToScreen import InfoToScreen
        builder = LearnerBuilder()
        builder = builder.buildRegressor()
        learner = builder.withTrainingDataFromARFFFile('refinement_strategy_%s'%self.example)\
        .withGrid().withLevel(3)\
        .withBorder(Types.BorderTypes.NONE)\
        .withSpecification().withIdentityOperator().withAdaptThreshold(0.001)\
        .withAdaptRate(1.0)\
        .withLambda(0.0001)\
        .withCGSolver().withImax(500)\
        .withStopPolicy().withAdaptiveItarationLimit(3)\
        .andGetResult()
        
        
        learner.specification.setBOperator(createOperationMultipleEval(learner.grid,
                    learner.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))
        
        while True: #repeat until policy says "stop"
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_STARTED)
            #learning step
            learner.alpha = learner.doLearningIteration(learner.dataContainer)
            learner.knowledge.update(learner.alpha)
#            a = learner.alpha.array().copy()
#            a.sort()
#            print a
            #compress grid
            if learner.iteration == 0:
                generator = learner.grid.createGridGenerator()
                functor = SurplusCoarseningFunctor(
                          learner.alpha,
                          generator.getNumberOfRemovablePoints(),
                          learner.specification.getAdaptThreshold())
                generator.coarsen(functor, learner.alpha)
            
            self.plotGrid(learner)
            
            storage = learner.grid.getStorage()
            grid_statistics = {}
            for i in xrange(storage.size()):
                point = storage.get(i)
#                print "[%d %d, %d %d]: %1.6f" % (point.getLevel(0),
#                      point.getIndex(0),point.getLevel(1),point.getIndex(1), 
#                      learner.alpha[i])
                key = (point.getLevel(0),point.getLevel(1))
                if grid_statistics.has_key(key):
                    grid_statistics[key] += 1
                else:
                    grid_statistics[key] = 1
            keys = sorted(grid_statistics.keys())
            
            l1max = max(keys, key=itemgetter(0))[0]
            l2max = max(keys, key=itemgetter(1))[1]
            m = zeros([l1max, l2max])
                
#            print learner.iteration,"Grid Statistics"
            for k in keys:
#                print "%s %s"%(k, grid_statistics[k])
                m[k[0]-1, k[1]-1] = grid_statistics[k]

            from pylab import pcolor, colorbar, savefig, clf, get_cmap
            pcolor(m,cmap=get_cmap('RdYlGn'))
            colorbar()
            savefig("%s%d_anova.png"%(self.example, learner.iteration))
            clf()
            
            del grid_statistics
            
            formatter = GridImageFormatter()
            formatter.serializeToFile(learner.grid, "%s%d_projections_anova.png"%(self.example, learner.iteration))
            clf()
            
            #calculate avg. error for training and test data and avg. for refine alpha
            learner.updateResults(learner.alpha, learner.dataContainer)
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            p_val = learner.trainAccuracy[-1] + learner.specification.getL()*np.sum(learner.alpha.array()**2)
            print "ANOVA iteration %d: %d grid points, %1.9f MSE, p* = %1.10f" % \
            (learner.iteration, storage.size(), learner.trainAccuracy[-1], p_val)           
            learner.iteration += 1
            if(learner.stopPolicy.isTrainingComplete(learner)): break
            
            #refine grid
            learner.notifyEventControllers(LearnerEvents.REFINING_GRID)
        
            pointsNum = learner.specification.getNumOfPointsToRefine( learner.grid.createGridGenerator().getNumberOfRefinablePoints() )
#            learner.grid.createGridGenerator().refine( SurplusRefinementFunctor(learner.errors, int(pointsNum), learner.specification.getAdaptThreshold()) )

            refiner = HashRefinement()
            functor = SurplusRefinementFunctor(learner.alpha, int(pointsNum), learner.specification.getAdaptThreshold())
            refinement_strategy = RefinementANOVAStrategy(functor)
            refiner.strategy_refine(learner.grid.getStorage(), refinement_strategy)
        formatter = GridFormatter()
        formatter.serializeToFile(learner.grid, "grid_anova.txt")
#        cg_solver = CGSolver()
#        x = cg_solver.solve()
#        print "Grid: ", learner.grid.serialize()



    def __testRefineStandard(self):
#        return
        from bin.pysgpp import HashRefinementBoundaries, RefinementANOVAStrategy, Grid, createOperationMultipleEval,\
        SurplusRefinementFunctor, HashRefinement, SurplusCoarseningFunctor
        from bin.learner.LearnerBuilder import LearnerBuilder
        from bin.learner import Types,  LearnerEvents
        from bin.data.DataContainer import DataContainer
#        from bin.controller.InfoToScreen import InfoToScreen
        builder = LearnerBuilder()
        builder = builder.buildRegressor()
        learner = builder.withTrainingDataFromARFFFile('refinement_strategy_%s'%self.example)\
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
#            a = learner.alpha.array().copy()
#            a.sort()
#            print a
            #compress grid
#            if learner.iteration == 0:
#                generator = learner.grid.createGridGenerator()
#                functor = SurplusCoarseningFunctor(
#                          learner.alpha,
#                          generator.getNumberOfRemovablePoints(),
#                          learner.specification.getAdaptThreshold())
#                generator.coarsen(functor, learner.alpha)
            
            storage = learner.grid.getStorage()
            grid_statistics = {}
            for i in xrange(storage.size()):
                point = storage.get(i)
#                print "[%d %d, %d %d]: %1.6f" % (point.getLevel(0),
#                      point.getIndex(0),point.getLevel(1),point.getIndex(1), 
#                      learner.alpha[i])
                key = (point.getLevel(0),point.getLevel(1))
                if grid_statistics.has_key(key):
                    grid_statistics[key] += 1
                else:
                    grid_statistics[key] = 1
            keys = sorted(grid_statistics.keys())
            
            l1max = max(keys, key=itemgetter(0))[0]
            l2max = max(keys, key=itemgetter(1))[1]
            m = zeros([l1max, l2max])
                
#            print learner.iteration,"Grid Statistics"
            for k in keys:
#                print "%s %s"%(k, grid_statistics[k])
                m[k[0]-1, k[1]-1] = grid_statistics[k]

            from pylab import pcolor, colorbar, savefig, clf, get_cmap
            pcolor(m,cmap=get_cmap('RdYlGn'))
            colorbar()
            savefig("%s%d_standard.png"%(self.example, learner.iteration))
            clf()
            
            del grid_statistics
            
            formatter = GridImageFormatter()
            formatter.serializeToFile(learner.grid, "%s%d_projections_standard.png"%(self.example, learner.iteration))
            clf()
            
            
            #calculate avg. error for training and test data and avg. for refine alpha
            learner.updateResults(learner.alpha, learner.dataContainer)
            learner.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            print "Standard iteration %d: %d grid points, %1.9f MSE" %(learner.iteration, storage.size(), learner.trainAccuracy[-1])
            learner.iteration += 1
            if(learner.stopPolicy.isTrainingComplete(learner)): break
            
            #refine grid
            learner.notifyEventControllers(LearnerEvents.REFINING_GRID)
        
            pointsNum = learner.specification.getNumOfPointsToRefine( learner.grid.createGridGenerator().getNumberOfRefinablePoints() )
            learner.grid.createGridGenerator().refine( SurplusRefinementFunctor(learner.alpha, int(pointsNum), learner.specification.getAdaptThreshold()) )
#            refiner = HashRefinement()
#            functor = SurplusRefinementFunctor(learner.errors, pointsNum, learner.specification.getAdaptThreshold())
#            refinement_strategy = RefinementANOVAStrategy(functor)
#            refiner.strategy_refine(learner.grid.getStorage(), refinement_strategy)
        formatter = GridFormatter()
        formatter.serializeToFile(learner.grid, "grid_normal.txt")
        
        
    def __testCG(self):
        return
        from bin.learner.LearnerBuilder import LearnerBuilder
        from bin.learner import Types,  LearnerEvents
        from bin.data.DataContainer import DataContainer
    #        from bin.controller.InfoToScreen import InfoToScreen
        builder = LearnerBuilder()
        builder = builder.buildRegressor()
        learner = builder.withTrainingDataFromARFFFile('refinement_strategy_%s'%self.example)\
        .withGrid().withLevel(3)\
        .withBorder(Types.BorderTypes.NONE)\
        .withSpecification().withIdentityOperator().withAdaptThreshold(0.001)\
        .withAdaptRate(1.0)\
        .withLambda(0.0000)\
        .withCGSolver().withImax(500)\
        .withStopPolicy().withAdaptiveItarationLimit(3)\
        .andGetResult()
        
        formatter = GridFormatter()
        grid_anova = formatter.deserializeFromFile("grid_anova.txt")
        grid_normal = formatter.deserializeFromFile("grid_normal.txt")
        
        grid = grid_anova
        
        linearSystem = DMSystemMatrix(grid,
                                       learner.dataContainer.getPoints(),
                                       learner.specification.getCOperator(),
                                       learner.specification.getL())
        size =  grid.getStorage().size() 
    
        alpha = DataVector(size)
        alpha.setAll( 0.0 )
        b1 = DataVector(size)
        linearSystem.generateb(learner.dataContainer.getValues(), b1)
        #calculates alphas
        learner.solver.solve(linearSystem, alpha, b1, learner.solver.getReuse(), 
                          True, learner.solver.getThreshold())
        
        
        matrix = np.zeros([size,size])
        eye_vector = DataVector(size)
        for i in xrange(size):
            col_vector = DataVector(size)
            
            eye_vector.setAll(0.0)
            eye_vector[i] = 1.0
            linearSystem.mult(eye_vector, col_vector)
            matrix[:,i] = col_vector.array()
        f = open("system_matrix_anova.pcl", 'w')
        pickle.dump(matrix, f)
        f.close()
        
        grid = grid_normal
        
        linearSystem = DMSystemMatrix(grid,
                                       learner.dataContainer.getPoints(),
                                       learner.specification.getCOperator(),
                                       learner.specification.getL())
        size =  grid.getStorage().size() 
    
        alpha = DataVector(size)
        alpha.setAll( 0.0 )
        b2 = DataVector(size)
        linearSystem.generateb(learner.dataContainer.getValues(), b2)
        #calculates alphas
        learner.solver.solve(linearSystem, alpha, b2, learner.solver.getReuse(), 
                          True, learner.solver.getThreshold())
        
        matrix = np.zeros([size,size])
        eye_vector = DataVector(size)
        for i in xrange(size):
            col_vector = DataVector(size)
            
            eye_vector.setAll(0.0)
            eye_vector[i] = 1.0
            linearSystem.mult(eye_vector, col_vector)
            matrix[:,i] = col_vector.array()
        f = open("system_matrix_normal.pcl", 'w')
        pickle.dump(matrix, f)
        f.close()
        
        print b1.toString()
        print b2.toString()
        
    
           
        
        
        
        
        
        
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()
