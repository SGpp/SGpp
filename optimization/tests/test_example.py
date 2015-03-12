# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import pysgpp
import objective_functions
import random
import math

class TestExample(unittest.TestCase):
    def setUp(self):
        """Initialize the test case."""
        # disable status output
        pysgpp.cvar.OptPrinterInstance.setVerbosity(-1)
        # disable multi-threading
        pysgpp.omp_set_num_threads(1)
    
    def testExample(self):
        """Test full example similar to c++_example.cpp."""
        random.seed(42)
        d = 2
        p = 3
        N = 100
        
        # test two simple objective functions
        fs = [objective_functions.ExampleFunction(), pysgpp.OptSphere(d)]
        f_gradients = [objective_functions.ExampleFunctionGradient(), objective_functions.SphereFunctionGradient(d)]
        f_hessians = [objective_functions.ExampleFunctionHessian(), objective_functions.SphereFunctionHessian(d)]
        # minima
        real_xopts = [[3.0/16.0 * math.pi, 3.0/14.0 * math.pi], [0.1, 0.1]]
        real_fopts = [-2.0, 0.0]
        # difference of global maximum/minimum
        function_ranges = [4.0, 81.0*d]
        
        grid = pysgpp.Grid.createModBsplineGrid(d, p)
        
        for f, f_gradient, f_hessian, real_xopt, real_fopt, function_range in \
                zip(fs, f_gradients, f_hessians, real_xopts, real_fopts, function_ranges):
            grid.getStorage().emptyStorage()
            gridgen = pysgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, 0.85)
            # generate grid
            self.assertTrue(gridgen.generate())
            
            # hierarchisation via OperationMultipleHierarchisation
            op_hier = pysgpp.createOperationMultipleHierarchisation(grid)
            alpha = pysgpp.DataVector(gridgen.getFunctionValues())
            op_hier.doHierarchisation(alpha)
            
            # create interpolant, gradient and Hessian
            ft = pysgpp.OptInterpolantFunction(d, grid, alpha)
            ft_gradient = pysgpp.OptInterpolantGradient(d, grid, alpha)
            ft_hessian = pysgpp.OptInterpolantHessian(d, grid, alpha)
            
            for i in range(100):
                # don't go near the boundary (should suffice)
                x = pysgpp.DataVector([random.uniform(0.2, 0.8) for t in range(d)])
                # test infinity norm of difference roughly
                self.assertLess(abs(f.eval(x) - ft.eval(x)), 0.25)
            
            # test all optimizers applied on function and interpolant
            optimizers = [pysgpp.OptGradientMethod(f, f_gradient),
                          pysgpp.OptNLCG(f, f_gradient),
                          pysgpp.OptNewton(f, f_hessian),
                          pysgpp.OptNelderMead(f),
                          pysgpp.OptRandomSearch(f),
                          pysgpp.OptDifferentialEvolution(f),
                          pysgpp.OptGradientMethod(ft, ft_gradient),
                          pysgpp.OptNLCG(ft, ft_gradient),
                          pysgpp.OptNewton(ft, ft_hessian),
                          pysgpp.OptNelderMead(ft),
                          pysgpp.OptRandomSearch(ft),
                          pysgpp.OptDifferentialEvolution(ft)]
            
            for optimizer in optimizers:
                xopt = pysgpp.DataVector(0)
                fopt = optimizer.optimize(xopt)
                self.assertEqual(len(xopt), d)
                # test distance of xopt in infinity norm
                for t in range(d): self.assertLessEqual(abs(xopt[t] - real_xopt[t]), 0.1)
                # test optimal function value
                self.assertAlmostEqual(fopt,
                        optimizer.getObjectiveFunction().eval(pysgpp.DataVector(xopt)))
                # allow 1% deviation of difference global maximum/minimum
                self.assertLessEqual(fopt - real_fopt, function_range / 100.0) 
