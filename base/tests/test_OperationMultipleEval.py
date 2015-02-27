# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/tbin/tlearner'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)


##
# Class with unittests for OperationMultipleEval
#
# @ingroup tests
#
class TestOperationMultipleEval(unittest.TestCase):

    ## 
    # Set up, create random DataVector and corresponding Python data structures.
    # @test DataVector::get(), DataVector::set()
    def setUp(self):
        pass

    ##
    # Constructors4.
    # @test DataVector::DataVector(size_t size), DataVector::DataVector(size_t size, size_t dim), DataVector::DataVector(DataVectorDefinition &DataVectorDef), DataVector::getSize(), DataVector::getDim(), DataVector::getSize()
    # @todo (pflueged) DataVector::DataVector(double *input, size_t size, size_t dim)
    def testConstructor(self):
        pass

    ##
    # Min, Max operations.
    # @test DataVector::min(int d), DataVector::max(int d), DataVector::minmax(int d, double *min, double *max), DataVector::min(), DataVector::max()
    def testOperationMultipleEval(self):
        from pysgpp import Grid
        from pysgpp import DataVector
        from pysgpp import DataMatrix
        from pysgpp import createOperationMultipleEval
        
        dim = 2
        level = 2
        grid = Grid.createLinearGrid(dim)
        grid.createGridGenerator().regular(level)
        gS = grid.getStorage()
        N = gS.size()
        alpha = DataVector([i + 1 for i in range(N)])
        
        points = [[0.5, 0.1], [0.3, 0.4], [0.9, 0.7]]
        
        numberDataPoints = 3
        
        result = DataVector([i for i in range(numberDataPoints)])
        
        dataset = DataMatrix(numberDataPoints, dim)
        for i in range(numberDataPoints):
            temp = DataVector(dim)
            for j in range(dim):
                print "p: ", points[i][j]
                temp[j] = points[i][j]
            dataset.setRow(i, temp)

        implementation = 7
        multiEvalOp = createOperationMultipleEval(grid, dataset, implementation)
        
        multiEvalOp.mult(alpha, result)        
        print 
        print "N: ", N
        print [alpha[i] for i in range(len(alpha))]
        print [result[i] for i in range(len(result))]
        #self.failUnlessEqual(quadOp.doQuadrature(alpha), qres)

# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()