# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


import unittest
import re
from pysgpp import DataVector, createOperationQuadrature


class TestQuadratureLinear(unittest.TestCase):
    ##
    # Test quadrature for a few test cases
    def testQuadrature(self):
        from pysgpp import Grid
        
        dim = 3
        level = 4
        grid = Grid.createLinearGrid(dim)
        grid.createGridGenerator().regular(level)
        gS = grid.getStorage()
        N = gS.size()
        alpha = DataVector([i for i in range(N)])

        # manual quadrature:
        qres = 0.0
        for i in range(N):
            lSum = gS.get(i).getLevelSum()
            qres += 2**(-lSum) * alpha[i]
        
        quadOp = createOperationQuadrature(grid)
        self.failUnlessEqual(quadOp.doQuadrature(alpha), qres)

class TestQuadratureMC(unittest.TestCase):
    ##
    # Test quadrature using Monte Carlo and compare against
    # sparse grid solution
    def testQuadratureMC(self):
        from pysgpp import Grid, OperationQuadratureMC
        
        dim = 3
        level = 3
        grid = Grid.createLinearGrid(dim)
        grid.createGridGenerator().regular(level)
        gS = grid.getStorage()
        N = gS.size()
        alpha = DataVector([i for i in range(N)])
        quadOp = createOperationQuadrature(grid)
        resDirect = quadOp.doQuadrature(alpha)

        # Monte Carlo quadrature
        opMC = OperationQuadratureMC(grid, 100000)
        resMC = opMC.doQuadrature(alpha)
        # only very coarse error check as randomized
        self.failUnlessAlmostEqual(resDirect, resMC, 0)

        
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()
