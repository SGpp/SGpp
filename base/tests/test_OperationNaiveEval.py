# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import pysgpp
import random

class TestOperationNaiveEval(unittest.TestCase):
    def testOperationNaiveEval(self):
        """Test naive evaluation operations (OperationNaiveEval)."""
        d = 2
        l = 4
        p = 3
        random.seed(42)
        
        # Test All The Grids!
        grids = [pysgpp.Grid.createBsplineGrid(d, p),
                 pysgpp.Grid.createBsplineTruncatedBoundaryGrid(d, p),
                 pysgpp.Grid.createBsplineClenshawCurtisGrid(d, p),
                 pysgpp.Grid.createModBsplineGrid(d, p),
                 pysgpp.Grid.createLinearGrid(d),
                 pysgpp.Grid.createLinearTruncatedBoundaryGrid(d),
                 pysgpp.Grid.createLinearClenshawCurtisGrid(d),
                 pysgpp.Grid.createModLinearGrid(d),
                 pysgpp.Grid.createWaveletGrid(d),
                 pysgpp.Grid.createWaveletTruncatedBoundaryGrid(d),
                 pysgpp.Grid.createModWaveletGrid(d)]
        
        bases = [pysgpp.SBsplineBase(p),
                 pysgpp.SBsplineBoundaryBase(p),
                 pysgpp.SBsplineClenshawCurtisBase(p),
                 pysgpp.SBsplineModifiedBase(p),
                 pysgpp.SLinearBase(),
                 pysgpp.SLinearBoundaryBase(),
                 pysgpp.SLinearClenshawCurtisBase(),
                 pysgpp.SLinearModifiedBase(),
                 pysgpp.SWaveletBase(),
                 pysgpp.SWaveletBoundaryBase(),
                 pysgpp.SWaveletModifiedBase()]
        
        for grid, basis in zip(grids, bases):
            # don't test gradients for linear function
            has_gradients = ("linear" not in grid.getType())
            
            # create regular sparse grid
            grid.createGridGenerator().regular(l)
            n = grid.getSize()
            alpha = pysgpp.DataVector(n)
            
            for i in range(n):
                gp = grid.getStorage().get(i)
                # don't forget to set the point distribution to Clenshaw-Curtis if necessary
                # (currently not done automatically)
                if grid.getType() in ["bsplineClenshawCurtis", "linearClenshawCurtis"]:
                    gp.setPointDistribution(pysgpp.HashGridIndex.ClenshawCurtis)
                x = pysgpp.DoubleVector([gp.getCoord(t) for t in range(d)])
                alpha[i] = random.gauss(0.0, 1.0)
            
            # create operations
            op = pysgpp.createOperationNaiveEval(grid)
            if has_gradients:
                op_dx = pysgpp.createOperationNaiveEvalGradient(grid)
                op_ddx = pysgpp.createOperationNaiveEvalHessian(grid)
                op_pdx = pysgpp.createOperationNaiveEvalPartialDerivative(grid)
            
            for k in range(20):
                # evaluate at random point
                x = pysgpp.DataVector([random.uniform(0.0, 1.0) for t in range(d)])
                fx = 0.0
                dfx = pysgpp.DoubleVector(d, 0.0)
                ddfx = pysgpp.DoubleVector(d*d, 0.0)
                
                for i in range(n):
                    # evaluate function by hand
                    gp = grid.getStorage().get(i)
                    val = alpha[i]
                    for t in range(d):
                        val *= basis.eval(gp.getLevel(t), gp.getIndex(t), x[t])
                    fx += val
                    
                    if not has_gradients: continue
                    
                    # evaluate gradient by hand
                    for j in range(d):
                        val = alpha[i]
                        for t in range(d):
                            if t == j:
                                val *= basis.evalDx(gp.getLevel(t), gp.getIndex(t), x[t])
                            else:
                                val *= basis.eval(gp.getLevel(t), gp.getIndex(t), x[t])
                        dfx[j] += val
                    
                    # evaluate Hessian by hand
                    for j in range(d):
                        for k in range(d):
                            val = alpha[i]
                            for t in range(d):
                                if (t == j) and (t == k):
                                    val *= basis.evalDxDx(gp.getLevel(t), gp.getIndex(t), x[t])
                                elif (t == j) or (t == k):
                                    val *= basis.evalDx(gp.getLevel(t), gp.getIndex(t), x[t])
                                else:
                                    val *= basis.eval(gp.getLevel(t), gp.getIndex(t), x[t])
                            ddfx[j*d+k] += val
                
                # test function evaluation
                fx2 = op.eval(alpha, x)
                self.assertAlmostEqual(fx, fx2, places=2)
                
                if has_gradients:
                    dfx2 = pysgpp.DataVector(d)
                    fx2 = op_dx.evalGradient(alpha, x, dfx2)
                    # test function evaluation
                    self.assertAlmostEqual(fx, fx2, places=2)
                    for t in range(d):
                        # test gradient evaluation
                        self.assertAlmostEqual(dfx[t], dfx2[t], places=2)
                        # test partial derivative evaluation
                        self.assertAlmostEqual(op_pdx.evalPartialDerivative(alpha, x, t), dfx[t])
                    
                    dfx2 = pysgpp.DataVector(d)
                    ddfx2 = pysgpp.DataMatrix(d, d)
                    fx2 = op_ddx.evalHessian(alpha, x, dfx2, ddfx2)
                    # test function evaluation
                    self.assertAlmostEqual(fx, fx2, places=2)
                    for t1 in range(d):
                        # test gradient evaluation
                        self.assertAlmostEqual(dfx[t1], dfx2[t1], places=2)
                        for t2 in range(d):
                            # test Hessian evaluation
                            self.assertAlmostEqual(ddfx[t1*d+t2], ddfx2.get(t1, t2), places=2) 
