###############################################################################
# Copyright (C) 2014 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
# @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

import math
import pysgpp
import random
import unittest

def testSLESystem(test_case, system, x, b):
    """Test sg::opt::sle::system::getMatrixEntry, isMatrixEntryNonZero and
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
    
    # A*x calculated by sg::opt
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

class TitleFunction(pysgpp.OptObjective):
    """Example objective function from the title of my Master's thesis."""
    def __init__(self):
        super(TitleFunction, self).__init__(2)
        # do not remove this, otherwise bad things will happen
        # (When exiting, Python tries to destruct clones made by clone(). But the clones
        # were already deleted by a SmartPointer in the C++ library, resulting in a double free
        # error. Somehow the reference circle created by the following line prevents Python from
        # doing the 2nd destruction, although it seems the problem has nothing to do with GC.)
        self.__this_prevents_destruction_by_python = self
    
    def eval(self, x):
        """Evaluates the function."""
        return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])
    
    def clone(self):
        """Clones the function object."""
        # create a new TitleFunction object and wrap it with a smart pointer
        # (The object is appended to an otherwise unneeded global list to prevent
        # Python from garbage collecting the object after returning from this function.)
        global sg_opt_clones
        sg_opt_clones.append(TitleFunction())
        return pysgpp.OptSmPtrObjective(sg_opt_clones[-1])

class TitleFunctionGradient(pysgpp.OptObjectiveGradient):
    """Gradient of TitleFunction."""
    def __init__(self):
        super(TitleFunctionGradient, self).__init__(2)
        # do not remove this, otherwise bad things will happen (see above)
        self.__this_prevents_destruction_by_python = self
    
    def evalGradient(self, x, gradient):
        """Evaluates the function gradient."""
        gradient[0] = 8.0 * math.cos(8.0 * x[0])
        gradient[1] = 7.0 * math.cos(7.0 * x[1])
        return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])
    
    def clone(self):
        """Clones the function object."""
        # see above
        global sg_opt_clones
        sg_opt_clones.append(TitleFunctionGradient())
        return pysgpp.OptSmPtrObjectiveGradient(sg_opt_clones[-1])

class TitleFunctionHessian(pysgpp.OptObjectiveHessian):
    """Gradient/Hessian of TitleFunction."""
    def __init__(self):
        super(TitleFunctionHessian, self).__init__(2)
        # do not remove this, otherwise bad things will happen (see above)
        self.__this_prevents_destruction_by_python = self
    
    def evalHessian(self, x, gradient, hessian):
        """Evaluates the function Hessian."""
        gradient[0] = 8.0 * math.cos(8.0 * x[0])
        gradient[1] = 7.0 * math.cos(7.0 * x[1])
        hessian.set(0, 0, -64.0 * math.sin(8.0 * x[0]))
        hessian.set(0, 1, 0.0)
        hessian.set(1, 0, 0.0)
        hessian.set(1, 1, -49.0 * math.sin(7.0 * x[1]))
        return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])
    
    def clone(self):
        """Clones the function object."""
        # see above
        global sg_opt_clones
        sg_opt_clones.append(TitleFunctionHessian())
        return pysgpp.OptSmPtrObjectiveHessian(sg_opt_clones[-1])

class SphereFunctionGradient(pysgpp.OptObjectiveGradient):
    def __init__(self, d):
        super(SphereFunctionGradient, self).__init__(d)
        # do not remove this, otherwise bad things will happen (see above)
        self.__this_prevents_destruction_by_python = self
    
    def evalGradient(self, x, gradient):
        """Evaluates the function Gradient."""
        d = self.getDimension()
        fx = 0.0
        for t in range(d):
            y = 10 * x[t] - 1
            fx += y*y
            gradient[t] = 20*y
        return fx
    
    def clone(self):
        """Clones the function object."""
        # see above
        global sg_opt_clones
        sg_opt_clones.append(SphereFunctionGradient(self.getDimension()))
        return pysgpp.OptSmPtrObjectiveGradient(sg_opt_clones[-1])

class SphereFunctionHessian(pysgpp.OptObjectiveHessian):
    def __init__(self, d):
        super(SphereFunctionHessian, self).__init__(d)
        # do not remove this, otherwise bad things will happen (see above)
        self.__this_prevents_destruction_by_python = self
    
    def evalHessian(self, x, gradient, hessian):
        """Evaluates the function Hessian."""
        d = self.getDimension()
        fx = 0.0
        for t in range(d):
            y = 10 * x[t] - 1
            fx += y*y
            gradient[t] = 20*y
            hessian.set(t, t, 200)
        return fx
    
    def clone(self):
        """Clones the function object."""
        # see above
        global sg_opt_clones
        sg_opt_clones.append(SphereFunctionHessian(self.getDimension()))
        return pysgpp.OptSmPtrObjectiveHessian(sg_opt_clones[-1])

class TestSGOpt(unittest.TestCase):
    """Test sg::opt classes and methods."""
    def setUp(self):
        """Initialize the test case."""
        # clear global clones list
        global sg_opt_clones
        sg_opt_clones = []
        # disable status output
        pysgpp.cvar.OptPrinterInstance.setVerbosity(-1)
        # disable multi-threading
        pysgpp.omp_set_num_threads(1)
    
    def testTestFunctions(self):
        """Test sg::opt::function::test::Test."""
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
            self.assertAlmostEqual(fopt, testfcn.eval(xopt))
            
            # test if xopt is minimal point for a sample of random points
            for i in range(1000):
                x = [random.uniform(0.0, 1.0) for t in range(d)]
                f = testfcn.eval(x)
                self.assertGreaterEqual(f, fopt)
    
    def testIterativeGridGenerators(self):
        """Test sg::opt::gridgen."""
        pysgpp.cvar.OptRNGInstance.setSeed(42)
        d = 2
        N = 200
        
        f = pysgpp.OptRosenbrock(d)
        f.generateDisplacement()
        
        # Test All The Grids!
        grids = [pysgpp.Grid.createBsplineGrid(d, 3),
                 pysgpp.Grid.createBsplineTrapezoidBoundaryGrid(d, 3),
                 pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3),
                 pysgpp.Grid.createModBsplineGrid(d, 3),
                 pysgpp.Grid.createLinearGrid(d),
                 pysgpp.Grid.createLinearTrapezoidBoundaryGrid(d),
                 pysgpp.Grid.createLinearClenshawCurtisGrid(d),
                 pysgpp.Grid.createModLinearGrid(d),
                 pysgpp.Grid.createWaveletGrid(d),
                 pysgpp.Grid.createWaveletTrapezoidBoundaryGrid(d),
                 pysgpp.Grid.createModWaveletGrid(d)]
        
        for grid in grids:
            # repeat for Ritter-Novak and linear surplus grid generation
            for k in range(2):
                # empty grid
                grid.getStorage().emptyStorage()
                if k == 0:
                    gridgen = pysgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, 0.85)
                else:
                    gridgen = pysgpp.OptIterativeGridGeneratorLinearSurplus(f, grid, N, 0.2)
                
                # generate grid
                self.assertTrue(gridgen.generate())
                
                # test grid size
                n = grid.getSize()
                self.assertLessEqual(n, N)
                
                # test size of function value vector
                function_values = gridgen.getFunctionValues()
                self.assertEqual(n, len(function_values))
                
                for i in range(n):
                    x = pysgpp.DoubleVector([grid.getStorage().get(i).abs(t) for t in range(d)])
                    # test function value
                    self.assertAlmostEqual(function_values[i], f.eval(x))
    
    def testOptimizers(self):
        """Test sg::opt::optimizer."""
        f = TitleFunction()
        f_gradient = TitleFunctionGradient()
        f_hessian = TitleFunctionHessian()
        
        # Test All The Optimizers!
        optimizers = [pysgpp.OptGradientMethod(f, f_gradient),
                      pysgpp.OptNLCG(f, f_gradient),
                      pysgpp.OptNewton(f, f_hessian),
                      pysgpp.OptNelderMead(f),
                      pysgpp.OptRandomSearch(f),
                      pysgpp.OptDifferentialEvolution(f)]
        
        for optimizer in optimizers:
            xopt = pysgpp.DoubleVector()
            # optimize
            fopt = optimizer.optimize(xopt)
            
            # test xopt and fopt
            self.assertEqual(len(xopt), 2)
            self.assertAlmostEqual(xopt[0], 3.0/16.0 * math.pi, places=2)
            self.assertAlmostEqual(xopt[1], 3.0/14.0 * math.pi, places=2)
            self.assertAlmostEqual(fopt, -2.0, places=2)
    
    def testSolvers(self):
        """Test sg::opt::sle::solver with sg::opt::sle::system::Full."""
        random.seed(42)
        # default solvers
        solvers = [pysgpp.OptBiCGStab(), pysgpp.OptAutoSolver()]
        
        # additional solvers if sg::opt was compiled with them
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
        """Test sg::opt::sle::system::Hierarchisation."""
        random.seed(42)
        d = 2
        l = 6
        
        f = TitleFunction()
        solver = pysgpp.OptAutoSolver()
        
        # Test All The Grids!
        grids = [pysgpp.Grid.createBsplineGrid(d, 3),
                 pysgpp.Grid.createBsplineTrapezoidBoundaryGrid(d, 3),
                 pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3),
                 pysgpp.Grid.createModBsplineGrid(d, 3),
                 pysgpp.Grid.createLinearGrid(d),
                 pysgpp.Grid.createLinearTrapezoidBoundaryGrid(d),
                 pysgpp.Grid.createLinearClenshawCurtisGrid(d),
                 pysgpp.Grid.createModLinearGrid(d),
                 pysgpp.Grid.createWaveletGrid(d),
                 pysgpp.Grid.createWaveletTrapezoidBoundaryGrid(d),
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
                if grid.getType() in ["BsplineClenshawCurtis", "linearClenshawCurtis"]:
                    gp.setPointDistribution(pysgpp.GridIndex.ClenshawCurtis)
                x = pysgpp.DoubleVector([gp.abs(t) for t in range(d)])
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
                self.assertLess(abs(f.eval(x) - ft.eval(x)), 0.1)
    
    def testNaiveEvalOperations(self):
        """Test naive evaluation operations created for sg::opt"""
        random.seed(42)
        d = 2
        l = 6
        p = 3
        
        f = TitleFunction()
        solver = pysgpp.OptAutoSolver()
        
        # Test All The Grids!
        grids = [pysgpp.Grid.createBsplineGrid(d, p),
                 pysgpp.Grid.createBsplineTrapezoidBoundaryGrid(d, p),
                 pysgpp.Grid.createBsplineClenshawCurtisGrid(d, p),
                 pysgpp.Grid.createModBsplineGrid(d, p),
                 pysgpp.Grid.createLinearGrid(d),
                 pysgpp.Grid.createLinearTrapezoidBoundaryGrid(d),
                 pysgpp.Grid.createLinearClenshawCurtisGrid(d),
                 pysgpp.Grid.createModLinearGrid(d),
                 pysgpp.Grid.createWaveletGrid(d),
                 pysgpp.Grid.createWaveletTrapezoidBoundaryGrid(d),
                 pysgpp.Grid.createModWaveletGrid(d)]
        
        bases = [pysgpp.SBsplineBase(p),
                 pysgpp.SBsplineBoundaryBase(p),
                 pysgpp.SBsplineClenshawCurtisBase(p),
                 pysgpp.SModBsplineBase(p),
                 pysgpp.SLinearBase(),
                 pysgpp.SLinearBoundaryBase(),
                 pysgpp.SLinearClenshawCurtisBase(),
                 pysgpp.SModLinearBase(),
                 pysgpp.SWaveletBase(),
                 pysgpp.SWaveletBoundaryBase(),
                 pysgpp.SModWaveletBase()]
        
        for grid, basis in zip(grids, bases):
            # don't test gradients for linear function
            has_gradients = ("linear" not in grid.getType())
            
            # create regular sparse grid
            grid.createGridGenerator().regular(l)
            n = grid.getSize()
            function_values = pysgpp.DoubleVector(n)
            
            for i in range(n):
                gp = grid.getStorage().get(i)
                # don't forget to set the point distribution to Clenshaw-Curtis if necessary
                # (currently not done automatically)
                if grid.getType() in ["BsplineClenshawCurtis", "linearClenshawCurtis"]:
                    gp.setPointDistribution(pysgpp.GridIndex.ClenshawCurtis)
                x = pysgpp.DoubleVector([gp.abs(t) for t in range(d)])
                function_values[i] = f.eval(x)
            
            # create hierarchisation system
            system = pysgpp.OptHierarchisationSystem(grid)
            alpha = pysgpp.DoubleVector(n)
            # solve system
            self.assertTrue(solver.solve(system, function_values, alpha))
            alpha_dv = pysgpp.DataVector(alpha)
            
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
                fx2 = op.eval(alpha_dv, x)
                self.assertAlmostEqual(fx, fx2, places=2)
                
                if has_gradients:
                    dfx2 = pysgpp.DataVector(d)
                    fx2 = op_dx.evalGradient(alpha_dv, x, dfx2)
                    # test function evaluation
                    self.assertAlmostEqual(fx, fx2, places=2)
                    for t in range(d):
                        # test gradient evaluation
                        self.assertAlmostEqual(dfx[t], dfx2[t], places=2)
                        # test partial derivative evaluation
                        self.assertAlmostEqual(op_pdx.evalPartialDerivative(alpha_dv, x, t),
                                               dfx[t])
                    
                    dfx2 = pysgpp.DataVector(d)
                    ddfx2 = pysgpp.DataMatrix(d, d)
                    fx2 = op_ddx.evalHessian(alpha_dv, x, dfx2, ddfx2)
                    # test function evaluation
                    self.assertAlmostEqual(fx, fx2, places=2)
                    for t1 in range(d):
                        # test gradient evaluation
                        self.assertAlmostEqual(dfx[t1], dfx2[t1], places=2)
                        for t2 in range(d):
                            # test Hessian evaluation
                            self.assertAlmostEqual(ddfx[t1*d+t2], ddfx2.get(t1, t2), places=2)
    
    def testExample(self):
        """Test full example similar to /sgpp/opt/example.cpp."""
        random.seed(42)
        d = 2
        p = 3
        N = 2000
        
        # test two simple objective functions
        fs = [TitleFunction(), pysgpp.OptSphere(d)]
        f_gradients = [TitleFunctionGradient(), SphereFunctionGradient(d)]
        f_hessians = [TitleFunctionHessian(), SphereFunctionHessian(d)]
        # minimas
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
            ft = pysgpp.OptInterpolant(d, grid, alpha)
            ft_gradient = pysgpp.OptInterpolantGradient(d, grid, alpha)
            ft_hessian = pysgpp.OptInterpolantHessian(d, grid, alpha)
            
            for i in range(100):
                # don't go near the boundary (should suffice)
                x = pysgpp.DoubleVector([random.uniform(0.2, 0.8) for t in range(d)])
                # test infinity norm of difference roughly
                self.assertLess(abs(f.eval(x) - ft.eval(x)), 0.1)
            
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
                xopt = pysgpp.DoubleVector()
                fopt = optimizer.optimize(xopt)
                self.assertEqual(len(xopt), d)
                # test distance of xopt in infinity norm
                for t in range(d): self.assertLessEqual(abs(xopt[t] - real_xopt[t]), 0.05)
                # test optimal function value
                self.assertAlmostEqual(fopt,
                        optimizer.getObjectiveFunction().eval(pysgpp.DoubleVector(xopt)))
                # allow 1% deviation of difference global maximum/minimum
                self.assertLessEqual(fopt - real_fopt, function_range / 100.0)
