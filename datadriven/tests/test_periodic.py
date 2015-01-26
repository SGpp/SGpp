# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

##############################################################################
# Copyright (C) 2014 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
# # @author Florian Zipperle (florian.zipperle@tum.de)


import unittest
import numpy as np
import math
from itertools import repeat, izip
             

class TestOperationEval(unittest.TestCase):
    def testEval1D(self):
        """Test the eval function for 1D grid"""
        from pysgpp import DataVector,DataMatrix, Grid
        
        grid = Grid.createPeriodicGrid(1)
        grid.createGridGenerator().regular(3)

        alpha = DataVector([1,] * grid.getStorage().size())
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([1])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0.5])))
        
        self.failUnlessEqual(2, grid.eval(alpha, DataVector([0.25])))
        self.failUnlessEqual(2, grid.eval(alpha, DataVector([0.75])))
        
        self.failUnlessEqual(2.5, grid.eval(alpha, DataVector([0.125])))
        
        alpha[0] = 0
        self.failUnlessEqual(0, grid.eval(alpha, DataVector([0])))
        self.failUnlessEqual(0, grid.eval(alpha, DataVector([1])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0.5])))
        
        self.failUnlessEqual(1.5, grid.eval(alpha, DataVector([0.25])))
        self.failUnlessEqual(1.5, grid.eval(alpha, DataVector([0.75])))
        
        self.failUnlessEqual(1.75, grid.eval(alpha, DataVector([0.125])))
        

    def testEval2D(self):
        """Test the eval function for 2D grid"""
        from pysgpp import DataVector,DataMatrix, Grid
        
        grid = Grid.createPeriodicGrid(2)
        grid.createGridGenerator().regular(3)

        alpha = DataVector([1,] * grid.getStorage().size())
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0,0])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([1,1])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0.5,0.5])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0,0.5])))
        
        self.failUnlessEqual(4, grid.eval(alpha, DataVector([0.25,0.25])))
        self.failUnlessEqual(4, grid.eval(alpha, DataVector([0.75,0.75])))
        
        
        self.failUnlessEqual(4.25, grid.eval(alpha, DataVector([0.125,0.125])))
        
        alpha[8] = 0
        self.failUnlessEqual(0, grid.eval(alpha, DataVector([0,0])))
        self.failUnlessEqual(0, grid.eval(alpha, DataVector([1,1])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0.5,0.5])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0,0.5])))
        
        self.failUnlessEqual(3.75, grid.eval(alpha, DataVector([0.25,0.25])))      
        
        self.failUnlessEqual(3.6875, grid.eval(alpha, DataVector([0.125,0.125])))
        
    def testEval3D(self):
        """Test the eval function for 3D grid"""
        from pysgpp import DataVector,DataMatrix, Grid
        
        grid = Grid.createPeriodicGrid(3)
        grid.createGridGenerator().regular(3)

        alpha = DataVector([1,] * grid.getStorage().size())
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0,0,0])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([1,1,1])))
        self.failUnlessEqual(1, grid.eval(alpha, DataVector([0.5,0.5,0.5])))
        self.failUnlessEqual(2, grid.eval(alpha, DataVector([0,0.5,0.25])))
        
        self.failUnlessEqual(7, grid.eval(alpha, DataVector([0.25,0.25,0.25])))        
        
        self.failUnlessEqual(6.25, grid.eval(alpha, DataVector([0.125,0.125,0.125])))
        
    def testMultEval1D(self):
        """Test the multiple eval function for 1D grid"""
        from pysgpp import DataVector,DataMatrix, Grid, OperationMultipleEval, createOperationMultipleEval
        
        grid = Grid.createPeriodicGrid(1)
        grid.createGridGenerator().regular(3)
        
        points = DataMatrix(5,1)
        points.set(0,0,0)
        points.set(1,0,1)
        points.set(2,0,0.5)
        points.set(3,0,0.25)
        points.set(4,0,0.125)
        
        eval = createOperationMultipleEval(grid,points)
        res = DataVector(5)
        
        alpha = DataVector([1,] * grid.getStorage().size())
        
        eval.mult(alpha,res)
        self.failUnlessEqual(1, res[0])
        self.failUnlessEqual(1, res[1])
        self.failUnlessEqual(1, res[2])       
        self.failUnlessEqual(2, res[3])     
        self.failUnlessEqual(2.5, res[4])
    
        alpha[0] = 0
        eval.mult(alpha,res)
        self.failUnlessEqual(0, res[0])
        self.failUnlessEqual(0, res[1])
        self.failUnlessEqual(1, res[2])       
        self.failUnlessEqual(1.5, res[3])     
        self.failUnlessEqual(1.75, res[4])
    
    def testMultEval2D(self):
        """Test the multiple eval function for 1D grid"""
        from pysgpp import DataVector,DataMatrix, Grid, OperationMultipleEval, createOperationMultipleEval
        
        grid = Grid.createPeriodicGrid(2)
        grid.createGridGenerator().regular(3)
        
        points = DataMatrix(3,2)
        points.set(0,0,0)
        points.set(0,1,0)
        
        points.set(1,0,0.5)
        points.set(1,1,0.5)
        
        points.set(2,0,0.125)
        points.set(2,1,0.125)
        
        eval = createOperationMultipleEval(grid,points)
        res = DataVector(3)
        
        alpha = DataVector([1,] * grid.getStorage().size())
        
        eval.mult(alpha,res)
        self.failUnlessEqual(1, res[0])
        self.failUnlessEqual(1, res[1])
        self.failUnlessEqual(4.25, res[2])       

        alpha[8] = 0
        eval.mult(alpha,res)
        self.failUnlessEqual(0, res[0])
        self.failUnlessEqual(1, res[1])
        self.failUnlessEqual(3.6875, res[2])  

class TestLTwoDotProduct(unittest.TestCase):
    def testLTwoExplicitD1L1(self):
        from pysgpp import OperationMatrixLTwoDotExplicitPeriodic, Grid, DataVector
        
        grid = Grid.createPeriodicGrid(1)
        grid.createGridGenerator().regular(1)
        
        ltwo = OperationMatrixLTwoDotExplicitPeriodic(grid)
        res = DataVector(2)
        
        alpha = DataVector([0,1])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.166666667, 0.333333333]), res)
        
        alpha = DataVector([1,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.333333333, 0.166666667]), res)
        
    def testLTwoExplicitD1L3(self):
        from pysgpp import OperationMatrixLTwoDotExplicitPeriodic, Grid, DataVector
        
        grid = Grid.createPeriodicGrid(1)
        grid.createGridGenerator().regular(3)
        
        ltwo = OperationMatrixLTwoDotExplicitPeriodic(grid)
        res = DataVector(8)
        
        alpha = DataVector([1,0,0,0,0,0,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.333333333333 , 0.166666666667 , 0.125 , 0.125 , 0.09375 , 0.03125 , 0.03125 , 0.09375]), res)
        
        alpha = DataVector([0,1,0,0,0,0,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.166666666667 , 0.333333333333 , 0.125 , 0.125 , 0.03125 , 0.09375 , 0.09375 , 0.03125]), res)
      
        alpha = DataVector([0,0,0,0,1,0,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.09375 , 0.03125 , 0.0625 , 0.0 , 0.0833333333333 , 0.0 , 0.0 , 0.0 ]), res)
        
        alpha = DataVector([0,0,0,0,0,0,0,1])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.09375 , 0.03125 , 0.0 , 0.0625 , 0.0 , 0.0 , 0.0 , 0.0833333333333]), res)
        
    def testLTwoExplicitD2L1(self):
        from pysgpp import OperationMatrixLTwoDotExplicitPeriodic, Grid, DataVector
        
        grid = Grid.createPeriodicGrid(2)
        grid.createGridGenerator().regular(1)
        
        ltwo = OperationMatrixLTwoDotExplicitPeriodic(grid)
        res = DataVector(4)
        
        alpha = DataVector([1,0,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.111111111111 , 0.0555555555556 , 0.0555555555556 , 0.0277777777778]), res)
        
        alpha = DataVector([0,1,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.0555555555556 , 0.111111111111 , 0.0277777777778 , 0.0555555555556]), res)
      
        alpha = DataVector([0,0,1,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.0555555555556 , 0.0277777777778 , 0.111111111111 , 0.0555555555556]), res)
        
        alpha = DataVector([0,0,0,1])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.0277777777778 , 0.0555555555556 , 0.0555555555556 , 0.111111111111]), res)
        
    def testLTwoExplicitD2L2(self):
        from pysgpp import OperationMatrixLTwoDotExplicitPeriodic, Grid, DataVector
        
        grid = Grid.createPeriodicGrid(2)
        grid.createGridGenerator().regular(2)
        
        ltwo = OperationMatrixLTwoDotExplicitPeriodic(grid)
        res = DataVector(12)
        
        alpha = DataVector([0,1,0,0,0,0,0,0,0,0,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.0555555555556 , 0.111111111111 , 0.0416666666667  , 0.0416666666667 ,0.0277777777778 ,  0.0208333333333 , 0.0208333333333 ,0.0555555555556 , 0.0416666666667, 0.0416666666667, 0.0208333333333 , 0.0208333333333  ]), res)
        
        
    def testLTwoD1L1(self):
        from pysgpp import OperationMatrixLTwoDotPeriodic, Grid, DataVector
        
        grid = Grid.createPeriodicGrid(1)
        grid.createGridGenerator().regular(1)
        
        ltwo = OperationMatrixLTwoDotPeriodic(grid.getStorage())
        res = DataVector(2)
        
        alpha = DataVector([0,1])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.166666667, 0.333333333]), res)
        
        alpha = DataVector([1,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.333333333, 0.166666667]), res)
        
    def testLTwoD1L3(self):
        from pysgpp import OperationMatrixLTwoDotPeriodic, Grid, DataVector
        
        grid = Grid.createPeriodicGrid(1)
        grid.createGridGenerator().regular(3)
        
        ltwo = OperationMatrixLTwoDotPeriodic(grid.getStorage())
        res = DataVector(8)
        
        alpha = DataVector([1,0,0,0,0,0,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.333333333333 , 0.166666666667 , 0.125 , 0.125 , 0.09375 , 0.03125 , 0.03125 , 0.09375]), res)
        
        alpha = DataVector([0,1,0,0,0,0,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.166666666667 , 0.333333333333 , 0.125 , 0.125 , 0.03125 , 0.09375 , 0.09375 , 0.03125]), res)
      
        alpha = DataVector([0,0,0,0,1,0,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.09375 , 0.03125 , 0.0625 , 0.0 , 0.0833333333333 , 0.0 , 0.0 , 0.0 ]), res)
        
        alpha = DataVector([0,0,0,0,0,0,0,1])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.09375 , 0.03125 , 0.0 , 0.0625 , 0.0 , 0.0 , 0.0 , 0.0833333333333]), res)
        
    def testLTwoD2L1(self):
        from pysgpp import OperationMatrixLTwoDotPeriodic, Grid, DataVector
        
        grid = Grid.createPeriodicGrid(2)
        grid.createGridGenerator().regular(1)
        
        ltwo = OperationMatrixLTwoDotPeriodic(grid.getStorage())
        res = DataVector(4)
        
        alpha = DataVector([1,0,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.111111111111 , 0.0555555555556 , 0.0555555555556 , 0.0277777777778]), res)
        
        alpha = DataVector([0,1,0,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.0555555555556 , 0.111111111111 , 0.0277777777778 , 0.0555555555556]), res)
      
        alpha = DataVector([0,0,1,0])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.0555555555556 , 0.0277777777778 , 0.111111111111 , 0.0555555555556]), res)
        
        alpha = DataVector([0,0,0,1])
        ltwo.mult(alpha, res)
        np.testing.assert_array_almost_equal(DataVector([0.0277777777778 , 0.0555555555556 , 0.0555555555556 , 0.111111111111]), res)
        
    def testLTwoD2L2(self):
        from pysgpp import OperationMatrixLTwoDotPeriodic, Grid, DataVector
        
        grid = Grid.createPeriodicGrid(2)
        grid.createGridGenerator().regular(2)
        
        ltwo = OperationMatrixLTwoDotPeriodic(grid.getStorage())
        res = DataVector(12)
        
        alpha = DataVector([0,1,0,0,0,0,0,0,0,0,0,0])
        ltwo.mult(alpha, res)
    
        np.testing.assert_array_almost_equal(DataVector([0.0555555555556 , 0.111111111111 , 0.0416666666667  , 0.0416666666667 ,0.0277777777778 ,  0.0208333333333 , 0.0208333333333 ,0.0555555555556 , 0.0416666666667, 0.0416666666667, 0.0208333333333 , 0.0208333333333  ]), res)
        
class TestLearnerDensityCluster(unittest.TestCase):
    def testExample1D_relative(self):
        from pysgpp import LearnerDensityCluster, DataMatrix, DataVector, Grid, SLESolverConfiguration, RegularGridConfiguration, Periodic, DensityBasedClusteringConfiguration
        data = [ 0.01, 0, 0.3, 0.7, 0.2, 0.78,  0.6, 0.82,  0.71, 0.18, 0.5]
        trainData = DataMatrix(10,1)
        for i in range(10):
            trainData.set(i, 0, data[i])
        
        gridConf = RegularGridConfiguration()
        gridConf.dim_ = 1;
        gridConf.level_ = 4;
        gridConf.type_ = Periodic;
    
        solvConf = SLESolverConfiguration()
        solvConf.eps_ = 0.0001;
        solvConf.maxIterations_ = 100;
        solvConf.threshold_ = -1.0;
        solvConf.type_ = 0 # sg::solver::CG;
    
        clustConf = DensityBasedClusteringConfiguration()
        clustConf.eps = 0.8;
        clustConf.numberOfNeighbors = 4;
        clustConf.thresholdType = 1;
    
        classes = DataVector(0)
        clust = LearnerDensityCluster()
        clust.setClusterConfiguration(clustConf)
        clust.train(trainData, classes, gridConf, solvConf, 0.01);   
        
        np.testing.assert_array_equal(clust.getClusterAssignments(), DataVector([0, 0, 1, 0, 0, 0, 2, 0, 0, 0]))
    
    def testExample1D_constant_classify(self):
        from pysgpp import LearnerDensityCluster, DataMatrix, DataVector, Grid, SLESolverConfiguration, RegularGridConfiguration, Periodic, DensityBasedClusteringConfiguration
        data = [ 0.01, 0, 0.3, 0.7, 0.2, 0.78,  0.6, 0.82,  0.71, 0.18, 0.5]
        trainData = DataMatrix(10,1)
        for i in range(10):
            trainData.set(i, 0, data[i])
        
        gridConf = RegularGridConfiguration()
        gridConf.dim_ = 1;
        gridConf.level_ = 4;
        gridConf.type_ = Periodic;
    
        solvConf = SLESolverConfiguration()
        solvConf.eps_ = 0.0001;
        solvConf.maxIterations_ = 100;
        solvConf.threshold_ = -1.0;
        solvConf.type_ = 0 # sg::solver::CG;
    
        clustConf = DensityBasedClusteringConfiguration()
        clustConf.eps = 10;
        clustConf.numberOfNeighbors = 4;
        clustConf.thresholdType = 0;
    
        classes = DataVector(0)
        clust = LearnerDensityCluster()
        clust.setClusterConfiguration(clustConf)
        clust.train(trainData, classes, gridConf, solvConf, 0.01);   
        res = clust.postprocessing(trainData,clust.getClusterAssignments(),2);
        
        np.testing.assert_array_equal(res, DataVector([ 0.,  0.,  1.,  0.,  1.,  0.,  0.,  0.,  0.,  1.]))
    
    def testExample1D_constant(self):
        from pysgpp import LearnerDensityCluster, DataMatrix, DataVector, Grid, SLESolverConfiguration, RegularGridConfiguration, Periodic, DensityBasedClusteringConfiguration
        data = [ 0.01, 0, 0.3, 0.7, 0.2, 0.78,  0.6, 0.82,  0.71, 0.18, 0.5]
        trainData = DataMatrix(10,1)
        for i in range(10):
            trainData.set(i, 0, data[i])
        
        gridConf = RegularGridConfiguration()
        gridConf.dim_ = 1;
        gridConf.level_ = 4;
        gridConf.type_ = Periodic;
    
        solvConf = SLESolverConfiguration()
        solvConf.eps_ = 0.0001;
        solvConf.maxIterations_ = 100;
        solvConf.threshold_ = -1.0;
        solvConf.type_ = 0 # sg::solver::CG;
    
        clustConf = DensityBasedClusteringConfiguration()
        clustConf.eps = 10;
        clustConf.numberOfNeighbors = 4;
        clustConf.thresholdType = 0;
    
        classes = DataVector(0)
        clust = LearnerDensityCluster()
        clust.setClusterConfiguration(clustConf)
        clust.train(trainData, classes, gridConf, solvConf, 0.01);   
        
        np.testing.assert_array_equal(clust.getClusterAssignments(), DataVector([0, 0, 1, 2, 1, 3, 4, 3, 2, 1]))   
        
    def testExample2D_minima(self):
        from pysgpp import LearnerDensityCluster, DataMatrix, DataVector, Grid, SLESolverConfiguration, RegularGridConfiguration, Periodic, DensityBasedClusteringConfiguration
        data = [ 0.5,0.5 , 0.2,0.2 , 0.2, 0.25 , 0.25,0.2, 0.8,0.85 ,0.8,0.8, 0.85,0.8]
        trainData = DataMatrix(7,2)
        for i in range(7):
            trainData.set(i, 0, data[2*i])
            trainData.set(i, 1, data[2*i+1])
        
        gridConf = RegularGridConfiguration()
        gridConf.dim_ = 2;
        gridConf.level_ = 4;
        gridConf.type_ = Periodic;
    
        solvConf = SLESolverConfiguration()
        solvConf.eps_ = 0.0001;
        solvConf.maxIterations_ = 100;
        solvConf.threshold_ = -1.0;
        solvConf.type_ = 0 # sg::solver::CG;
    
        clustConf = DensityBasedClusteringConfiguration()
        clustConf.eps = 0.1;
        clustConf.numberOfNeighbors = 3;
        clustConf.thresholdType = 3;
        clustConf.threshold = 0.3
    
        classes = DataVector(0)
        clust = LearnerDensityCluster()
        clust.setClusterConfiguration(clustConf)
        clust.train(trainData, classes, gridConf, solvConf, 0.01);   
        
        np.testing.assert_array_equal(clust.getClusterAssignments(), DataVector([ 0.,  0.,  0.,  0.,  1.,  1.,  1.]))           
if __name__ == "__main__":
    unittest.main()   