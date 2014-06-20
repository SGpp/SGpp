##############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
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

                
if __name__ == "__main__":
    unittest.main()   
