# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import pysgpp
import random

class TestTestFunctions(unittest.TestCase):
    def setUp(self):
        """Initialize the test case."""
        # disable status output
        pysgpp.cvar.OptPrinterInstance.setVerbosity(-1)
        # disable multi-threading
        pysgpp.omp_set_num_threads(1)
    
    def testTestFunctions(self):
        """Test SGPP::optimization::test_functions::TestFunction."""
        random.seed(42)
        pysgpp.cvar.OptRNGInstance.setSeed(42)
        
        d = 10
        # Test All The Testfunctions!
        testfcns = [pysgpp.OptAckley(d),
                    pysgpp.OptBeale(),
                    pysgpp.OptBranin(),
                    pysgpp.OptEasom(),
                    pysgpp.OptEggholder(),
                    pysgpp.OptGoldsteinPrice(),
                    pysgpp.OptGriewank(d),
                    pysgpp.OptHartman3(),
                    pysgpp.OptHartman6(),
                    pysgpp.OptHimmelblau(),
                    pysgpp.OptHoelderTable(),
                    pysgpp.OptMichalewicz(),
                    pysgpp.OptMladineo(),
                    pysgpp.OptRastrigin(d),
                    pysgpp.OptRosenbrock(d),
                    pysgpp.OptSchwefel(d),
                    pysgpp.OptSHCB(),
                    pysgpp.OptSphere(d)]
        
        for testfcn in testfcns:
            d = testfcn.getDimension()
            # displace function randomly
            testfcn.generateDisplacement()
            
            # generate random point, displace and reverse
            xl1 = [random.uniform(0.0, 1.0) for t in range(d)]
            x = pysgpp.DoubleVector(xl1)
            testfcn.displaceVector(x)
            f1 = testfcn.evalUndisplaced(x)
            testfcn.reverseDisplaceVector(x)
            f2 = testfcn.eval(x)
            xl2 = list(x)
            # test displaceVector/reverseDisplaceVector
            for t in range(d): self.assertAlmostEqual(xl1[t], xl2[t])
            # test eval/evalUndisplaced
            self.assertAlmostEqual(f1, f2)
            
            xopt = pysgpp.DoubleVector()
            fopt = testfcn.getOptimalPoint(xopt)
            # test minimal point
            self.assertEqual(len(xopt), d)
            for t in range(d):
                self.assertGreaterEqual(xopt[t], 0.0)
                self.assertLessEqual(xopt[t], 1.0)
            
            if pysgpp.cvar.USING_DOUBLE_PRECISION:
                self.assertAlmostEqual(fopt, testfcn.eval(xopt))
            else:
                self.assertAlmostEqual(fopt, testfcn.eval(xopt), places=6)
            
            # test if xopt is minimal point for a sample of random points
            for i in range(1000):
                x = [random.uniform(0.0, 1.0) for t in range(d)]
                f = testfcn.eval(x)
                self.assertGreaterEqual(f, fopt)
