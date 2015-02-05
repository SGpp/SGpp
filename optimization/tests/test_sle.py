# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import pysgpp
import objective_functions
import random
import math

def testSLESystem(test_case, system, x, b):
    """Test SGPP::optimization::sle::system::getMatrixEntry, isMatrixEntryNonZero and
    matrixVectorMultiplication. Returns system matrix as pysgpp.DataMatrix.
    """
    n = len(x)
    A = pysgpp.DataMatrix(n, n)
    # A*x calculated directly
    Ax = pysgpp.DoubleVector(n, 0.0)
    for i in range(n):
        for j in range(n):
            Aij = system.getMatrixEntry(i, j)
            A.set(i, j, Aij)
            Ax[i] += Aij * x[j]
            # test isMatrixEntryNonZero
            test_case.assertEqual(system.isMatrixEntryNonZero(i, j), Aij != 0)
    
    # A*x calculated by SGPP::opt
    Ax2 = pysgpp.DoubleVector()
    system.matrixVectorMultiplication(x, Ax2)
    for i in range(n): test_case.assertAlmostEqual(Ax[i], Ax2[i])
    return A

def testSLESolution(test_case, A, x, b):
    """Test A*x = b."""
    n = len(x)
    r_norm_squared = 0.0
    b_norm_squared = 0.0
    for i in range(n):
        ri = b[i]
        for j in range(n): ri -= A.get(i, j) * x[j]
        r_norm_squared += ri*ri
        b_norm_squared += b[i]*b[i]
    # test relative residual
    test_case.assertLess(math.sqrt(r_norm_squared / b_norm_squared), 1e-6)

class TestSLE(unittest.TestCase):
    def setUp(self):
        """Initialize the test case."""
        # disable status output
        pysgpp.cvar.OptPrinterInstance.setVerbosity(-1)
        # disable multi-threading
        pysgpp.omp_set_num_threads(1)
    
    def testSolvers(self):
        """Test SGPP::optimization::sle::solver with SGPP::opt::sle::system::Full."""
        random.seed(42)
        # default solvers
        solvers = [pysgpp.OptBiCGStab(),
                   pysgpp.OptGaussianElimination(),
                   pysgpp.OptAutoSolver()]
        
        # additional solvers if SGPP::opt was compiled with them
        if pysgpp.cvar.ARMADILLO_ENABLED:
            solvers.append(pysgpp.OptArmadillo())
        if pysgpp.cvar.EIGEN_ENABLED:
            solvers.append(pysgpp.OptEigen())
        if pysgpp.cvar.GMMPP_ENABLED:
            solvers.append(pysgpp.OptGmmpp())
        if pysgpp.cvar.UMFPACK_ENABLED:
            solvers.append(pysgpp.OptUMFPACK())
        
        # test various SLE dimensions
        for n in (list(range(1, 11)) + [20, 50, 100, 200]):
            # generate random matrix and RHS
            A = pysgpp.DataMatrix(n, n)
            b = pysgpp.DoubleVector(n)
            
            for i in range(n):
                b[i] = random.uniform(-1.0, 1.0)
                for j in range(n):
                    A.set(i, j, random.uniform(-1.0, 1.0))
            
            # full SLE
            system = pysgpp.OptFullSystem(A)
            
            #print "n = " + str(n)
            #print "A = [",
            #for i in range(n):
            #    if i > 0: print ";"
            #    for j in range(n):
            #        if j > 0: print ", ",
            #        print A.get(i, j),
            #print "]"
            #
            #print "b = [",
            #for i in range(n):
            #    if i > 0: print "; ",
            #    print b[i],
            #print "]"
            
            for solver in solvers:
                if (type(solver) is pysgpp.OptBiCGStab) and (n > 20):
                    # BiCGStab is really weak and can't solve bigger systems
                    # (maybe there's a bug in there, but it should only be used for Newton's
                    # method of optimization - for hierarchisation, only external solvers should
                    # be used)
                    continue
                
                x = pysgpp.DoubleVector(n)
                # solve system
                self.assertTrue(solver.solve(system, b, x))
                # test solution
                testSLESolution(self, A, x, b)
            
            # test system with last obtained solution (should suffice)
            testSLESystem(self, system, x, b)
    
    def testHierarchisation(self):
        """Test SGPP::optimization::sle::system::Hierarchisation."""
        random.seed(42)
        d = 2
        l = 4
        
        f = objective_functions.TitleFunction()
        solver = pysgpp.OptAutoSolver()
        
        # Test All The Grids!
        grids = [pysgpp.Grid.createBsplineGrid(d, 3),
                 pysgpp.Grid.createBsplineTruncatedBoundaryGrid(d, 3),
                 pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3),
                 pysgpp.Grid.createModBsplineGrid(d, 3),
                 pysgpp.Grid.createLinearGrid(d),
                 pysgpp.Grid.createLinearTruncatedBoundaryGrid(d),
                 pysgpp.Grid.createLinearClenshawCurtisGrid(d),
                 pysgpp.Grid.createModLinearGrid(d),
                 pysgpp.Grid.createWaveletGrid(d),
                 pysgpp.Grid.createWaveletTruncatedBoundaryGrid(d),
                 pysgpp.Grid.createModWaveletGrid(d)]
        
        for grid in grids:
            # generate regular sparse grid
            grid.createGridGenerator().regular(l)
            n = grid.getSize()
            function_values = pysgpp.DoubleVector(n)
            
            for i in range(n):
                gp = grid.getStorage().get(i)
                # don't forget to set the point distribution to Clenshaw-Curtis if necessary
                # (currently not done automatically)
                if grid.getType() in ["bsplineClenshawCurtis", "linearClenshawCurtis"]:
                    gp.setPointDistribution(pysgpp.HashGridIndex.ClenshawCurtis)
                x = pysgpp.DoubleVector([gp.getCoord(t) for t in range(d)])
                function_values[i] = f.eval(x)
            
            # create hierarchisation system
            system = pysgpp.OptHierarchisationSystem(grid)
            alpha = pysgpp.DoubleVector(n)
            # solve system
            self.assertTrue(solver.solve(system, function_values, alpha))
            
            # test system
            A = testSLESystem(self, system, alpha, function_values)
            # test solution
            testSLESolution(self, A, alpha, function_values)
            
            # create interpolant
            ft = pysgpp.OptInterpolant(d, grid, pysgpp.DataVector(alpha))
            for i in range(100):
                # don't go near the boundary (should suffice)
                x = pysgpp.DoubleVector([random.uniform(0.2, 0.8) for t in range(d)])
                # test infinity norm of difference roughly
                self.assertLess(abs(f.eval(x) - ft.eval(x)), 0.3)
