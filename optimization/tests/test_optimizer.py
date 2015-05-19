# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import pysgpp
import objective_functions as objfcns
import random
import math

class TestOptimizer(unittest.TestCase):
    def setUp(self):
        """Initialize the test case."""
        # disable status output
        pysgpp.cvar.OptPrinterInstance.setVerbosity(-1)
        # disable multi-threading
        pysgpp.omp_set_num_threads(1)
    
    def testUnconstrainedOptimizers(self):
        """Test unconstrained optimizers in SGPP::optimization::optimizer."""
        f = objfcns.ExampleFunction()
        f_gradient = objfcns.ExampleFunctionGradient()
        f_hessian = objfcns.ExampleFunctionHessian()
        
        # Test All The Optimizers!
        N = 1000
        optimizers = [pysgpp.OptGradientDescent(f, f_gradient, N),
                      pysgpp.OptNLCG(f, f_gradient, N),
                      pysgpp.OptNewton(f, f_hessian, N),
                      pysgpp.OptAdaptiveGradientDescent(f, f_gradient, N),
                      pysgpp.OptAdaptiveNewton(f, f_hessian, N),
                      pysgpp.OptBFGS(f, f_gradient, N),
                      pysgpp.OptRprop(f, f_gradient, N),
                      pysgpp.OptNelderMead(f, N),
                      pysgpp.OptMultiStart(f, N),
                      pysgpp.OptDifferentialEvolution(f, N)]
        
        for optimizer in optimizers:
            xopt = pysgpp.DataVector(0)
            # set starting point
            x0 = pysgpp.DataVector(2)
            x0[0] = 0.8
            x0[1] = 0.5
            optimizer.setStartingPoint(x0)
            # optimize
            fopt = optimizer.optimize(xopt)
            print optimizer
            print xopt
            print fopt
            
            # test xopt and fopt
            self.assertEqual(len(xopt), 2)
            self.assertAlmostEqual(xopt[0], 3.0/16.0 * math.pi, places=3)
            self.assertAlmostEqual(xopt[1], 3.0/14.0 * math.pi, places=3)
            self.assertAlmostEqual(fopt, -2.0, places=6)
    
    def testConstrainedOptimizers(self):
        """Test constrained optimizers in SGPP::optimization::optimizer."""
        N = 10000
        
        for i in range(3):
            if i == 0:
                d = 5
                xopt_real = d * [1.0 / math.sqrt(d)]
                fopt_real = 1.0
                x0 = pysgpp.DataVector(d, 0.5)
                f = objfcns.G3ObjectiveFunction(d)
                f_gradient = objfcns.G3ObjectiveGradient(d)
                g = pysgpp.cvar.OptEmptyConstraintFunctionInstance
                g_gradient = pysgpp.cvar.OptEmptyConstraintGradientInstance
                h = objfcns.G3ConstraintFunction(d)
                h_gradient = objfcns.G3ConstraintGradient(d)
                optimizers = [pysgpp.OptSquaredPenalty(f, f_gradient,
                                                       g, g_gradient,
                                                       h, h_gradient, N),
                              #pysgpp.OptLogBarrier(f, f_gradient,
                              #                     g, g_gradient, N),
                              pysgpp.OptAugmentedLagrangian(f, f_gradient,
                                                            g, g_gradient,
                                                            h, h_gradient, N)]
            elif i == 1:
                d = 2
                xopt_real = [(289.0/19.0 - 13.0) / (100.0 - 13.0),
                             ((1000.0 - math.sqrt(691239.0)) / 200.0 - 0.0) /
                             (100.0 -  0.0)]
                fopt_real = -(3000.0 + math.sqrt(691239.0))**3 / 8000000.0 + \
                            549353259.0 / 8000000.0
                x0 = pysgpp.DataVector(2, 0.5)
                f = objfcns.G6ObjectiveFunction()
                f_gradient = objfcns.G6ObjectiveGradient()
                g = objfcns.G6ConstraintFunction()
                g_gradient = objfcns.G6ConstraintGradient()
                h = pysgpp.cvar.OptEmptyConstraintFunctionInstance
                h_gradient = pysgpp.cvar.OptEmptyConstraintGradientInstance
                optimizers = [pysgpp.OptSquaredPenalty(f, f_gradient,
                                                       g, g_gradient,
                                                       h, h_gradient, N),
                              pysgpp.OptLogBarrier(f, f_gradient,
                                                   g, g_gradient, N),
                              pysgpp.OptAugmentedLagrangian(f, f_gradient,
                                                            g, g_gradient,
                                                            h, h_gradient, N)]
            else:
                d = 2
                xopt_real = [1.2279713 / 10.0, 4.2453733 / 10.0]
                fopt_real = 0.095825
                x0 = pysgpp.DataVector(2, 0.5)
                f = objfcns.G8ObjectiveFunction()
                f_gradient = objfcns.G8ObjectiveGradient()
                g = objfcns.G8ConstraintFunction()
                g_gradient = objfcns.G8ConstraintGradient()
                h = pysgpp.cvar.OptEmptyConstraintFunctionInstance
                h_gradient = pysgpp.cvar.OptEmptyConstraintGradientInstance
                optimizers = [pysgpp.OptSquaredPenalty(f, f_gradient,
                                                       g, g_gradient,
                                                       h, h_gradient, N),
                              pysgpp.OptLogBarrier(f, f_gradient,
                                                   g, g_gradient, N),
                              pysgpp.OptAugmentedLagrangian(f, f_gradient,
                                                            g, g_gradient,
                                                            h, h_gradient, N)]
            
            for optimizer in optimizers:
                print f
                print optimizer
                xopt = pysgpp.DataVector(0)
                # set starting point
                optimizer.setStartingPoint(x0)
                # optimize
                fopt = optimizer.optimize(xopt)
                print xopt
                print fopt
                
                # test xopt and fopt
                #self.assertEqual(len(xopt), d)
                #for t in range(d):
                #    self.assertAlmostEqual(xopt[t], xopt_real[t], places=3)
                #self.assertAlmostEqual(fopt, fopt_real, places=6)
