###############################################################################
# Copyright (C) 2014 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
# @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

import math
import pysgpp
import unittest

def testUniformUnmodifiedLinearFunctions(test_case, basis):
    """Test unmodified linear basis functions for level >= 1."""
    test_case.assertEqual(basis.eval(1, 1, 0.5),   1.0)
    test_case.assertEqual(basis.eval(1, 1, 0.75),  0.5)
    test_case.assertEqual(basis.eval(1, 1, 0.875), 0.25)
    test_case.assertEqual(basis.eval(1, 1, 0.0),   0.0)
    test_case.assertEqual(basis.eval(1, 1, 1.0),   0.0)
    test_case.assertEqual(basis.eval(2, 1, 0.75),  0.0)
    test_case.assertEqual(basis.eval(2, 3, 0.75),  1.0)
    test_case.assertEqual(basis.eval(3, 1, 0.0),   0.0)
    test_case.assertEqual(basis.eval(3, 1, 0.125), 1.0)
    test_case.assertEqual(basis.eval(3, 1, 0.25),  0.0)

def testLinearLevelZeroFunctions(test_case, basis):
    """Test boundary linear basis functions (level 0)."""
    for i in range(2):
        for x in [j/4.0 for j in range(5)]:
            test_case.assertEqual(basis.eval(0, i, x), i*x + (1-i)*(1.0-x))

def testModLinearFunctions(test_case, basis):
    """Test modified linear basis functions."""
    test_case.assertEqual(basis.eval(1, 1, 0.5),   1.0)
    test_case.assertEqual(basis.eval(1, 1, 0.75),  1.0)
    test_case.assertEqual(basis.eval(1, 1, 0.875), 1.0)
    test_case.assertEqual(basis.eval(1, 1, 0.0),   1.0)
    test_case.assertEqual(basis.eval(1, 1, 1.0),   1.0)
    test_case.assertEqual(basis.eval(2, 1, 0.0),   2.0)
    test_case.assertEqual(basis.eval(2, 1, 0.75),  0.0)
    test_case.assertEqual(basis.eval(2, 3, 0.75),  1.0)
    test_case.assertEqual(basis.eval(2, 3, 1.0),   2.0)
    test_case.assertEqual(basis.eval(3, 1, 0.0),   2.0)
    test_case.assertEqual(basis.eval(3, 1, 0.125), 1.0)
    test_case.assertEqual(basis.eval(3, 1, 0.25),  0.0)

def ccKnot(l, i):
    """Return Clenshaw-Curtis knot with given level and index."""
    return 0.5 * (math.cos(math.pi * (1.0 - i / 2.0**l)) + 1.0)

def testLinearClenshawCurtisFunctions(test_case, basis):
    """Test Clenshaw-Curtis linear basis functions for level >= 1."""
    test_case.assertEqual(basis.eval(1, 1, 0.5),   1.0)
    test_case.assertEqual(basis.eval(1, 1, 0.75),  0.5)
    test_case.assertEqual(basis.eval(1, 1, 0.875), 0.25)
    test_case.assertEqual(basis.eval(1, 1, 0.0),   0.0)
    test_case.assertEqual(basis.eval(1, 1, 1.0),   0.0)
    test_case.assertAlmostEqual(basis.eval(2, 1, 0.75),         0.0)
    test_case.assertAlmostEqual(basis.eval(2, 3, 0.75),         0.25 / (ccKnot(2, 3) - 0.5))
    test_case.assertAlmostEqual(basis.eval(2, 3, ccKnot(2, 3)), 1.0)
    test_case.assertAlmostEqual(basis.eval(3, 1, 0.0),          0.0)
    test_case.assertAlmostEqual(basis.eval(3, 1, 0.125),
            1.0 - (0.125 - ccKnot(3, 1)) / (ccKnot(3, 2) - ccKnot(3, 1)))
    test_case.assertAlmostEqual(basis.eval(3, 1, ccKnot(3, 1)), 1.0)
    test_case.assertAlmostEqual(basis.eval(3, 1, 0.25),         0.0)

def testDerivatives(test_case, basis, deg=2, start_level=1, max_discontinuities_count=0):
    """Test derivatives (up to order deg, max. 2) of basis functions for level >= start_level.
    Allow for max_discontinuities_count discontinuities (e.g. for Wavelets which are cut off).
    """
    dx = 1e-10
    n = 2
    
    # levels
    for l in range(start_level, 4):
        # indices
        for i in range(1, 2**l, 2):
            f = lambda y: basis.eval(l, i, y)
            
            if deg >= 1:
                # test first derivative at boundary (central difference quotient)
                df = lambda y: basis.evalDx(l, i, y)
                test_case.assertAlmostEqual((f(2*dx) - f(0.0)) / (2*dx), df(dx), places=n)
                test_case.assertAlmostEqual((f(1.0) - f(1.0-2*dx)) / (2*dx), df(1.0-dx), places=n)
            if deg >= 2:
                # test second derivative at boundary (central difference quotient)
                ddf = lambda y: basis.evalDxDx(l, i, y)
                test_case.assertAlmostEqual((df(2*dx) - df(0.0)) / (2*dx), ddf(dx), places=n)
                test_case.assertAlmostEqual((df(1.0) - df(1.0-2*dx)) / (2*dx), ddf(1.0-dx),
                                            places=n)
            
            # count discontinuities
            discontinuities = 0
            for x in [j/100.0 for j in range(1, 100)]:
                if abs(f(x+dx) - f(x-dx)) > 1e5 * dx:
                    # discontinuity found
                    discontinuities += 1
                else:
                    # test derivatives only if the function is continuous
                    if deg >= 1:
                        # test first derivative (central difference quotient)
                        test_case.assertAlmostEqual((f(x+dx) - f(x-dx)) / (2*dx), df(x),
                                                    places=n)
                    if deg >= 2:
                        # test second derivative (central difference quotient)
                        test_case.assertAlmostEqual((df(x+dx) - df(x-dx)) / (2*dx), ddf(x),
                                                    places=n)
            test_case.assertLessEqual(discontinuities, max_discontinuities_count)

def testBsplineProperties(test_case, basis, start_level=1):
    """Test basic B-spline properties (mixed monotonicity, bounds) for level >= start_level."""
    for l in range(start_level, 4):
        for i in range(1, 2**l, 2):
            # rising at the beginning
            falling = False
            fx = None
            for x in [j/100.0 for j in range(101)]:
                fx_new = basis.eval(l, i, x)
                if fx is not None:
                    if falling:
                        # hope we're still falling
                        test_case.assertLessEqual(fx_new, fx)
                    elif fx_new < fx:
                        # we're now falling (and weren't until now)
                        falling = True
                fx = fx_new
                
                # test bounds
                test_case.assertGreaterEqual(fx, 0.0)
                test_case.assertLessEqual(fx, 1.0)

class TestBasis(unittest.TestCase):
    """Test those basis functions which are supported by sg::opt."""
    def testLinearBasis(self):
        """Test linear Noboundary basis."""
        basis = pysgpp.SLinearBase()
        testUniformUnmodifiedLinearFunctions(self, basis)
        testDerivatives(self, basis, deg=0)
    
    def testLinearBoundaryBasis(self):
        """Test linear Boundary basis."""
        basis = pysgpp.SLinearBoundaryBase()
        testLinearLevelZeroFunctions(self, basis)
        testUniformUnmodifiedLinearFunctions(self, basis)
        testDerivatives(self, basis, deg=0, start_level=0)
    
    def testLinearClenshawCurtisBasis(self):
        """Test linear Clenshaw-Curtis basis."""
        basis = pysgpp.SLinearClenshawCurtisBase()
        testLinearLevelZeroFunctions(self, basis)
        testLinearClenshawCurtisFunctions(self, basis)
        testDerivatives(self, basis, deg=0, start_level=0)
    
    def testModLinearBasis(self):
        """Test modified linear basis."""
        basis = pysgpp.SModLinearBase()
        testModLinearFunctions(self, basis)
        testDerivatives(self, basis, deg=0)
    
    def testBsplineBasis(self):
        """Test B-spline Noboundary basis."""
        basis = pysgpp.SBsplineBase(1)
        testUniformUnmodifiedLinearFunctions(self, basis)
        
        for p in range(6):
            basis = pysgpp.SBsplineBase(p)
            testBsplineProperties(self, basis)
            testDerivatives(self, basis, deg=max(basis.getDegree() - 1, 0))
    
    def testBsplineBoundaryBasis(self):
        """Test B-spline Boundary basis."""
        basis = pysgpp.SBsplineBoundaryBase(1)
        testLinearLevelZeroFunctions(self, basis)
        testUniformUnmodifiedLinearFunctions(self, basis)
        
        for p in range(6):
            basis = pysgpp.SBsplineBase(p)
            testBsplineProperties(self, basis, 0)
            testDerivatives(self, basis, deg=max(basis.getDegree() - 1, 0), start_level=0)
    
    def testBsplineClenshawCurtisBasis(self):
        """Test B-spline Clenshaw-Curtis basis."""
        basis = pysgpp.SBsplineClenshawCurtisBase(1)
        testLinearLevelZeroFunctions(self, basis)
        testLinearClenshawCurtisFunctions(self, basis)
        
        for p in range(6):
            basis = pysgpp.SBsplineBase(p)
            testBsplineProperties(self, basis, 0)
            testDerivatives(self, basis, deg=max(basis.getDegree() - 1, 0), start_level=0)
    
    def testModBsplineBasis(self):
        """Test modified B-spline basis."""
        basis = pysgpp.SModBsplineBase(1)
        testModLinearFunctions(self, basis)
        
        for p in range(6):
            basis = pysgpp.SBsplineBase(p)
            testBsplineProperties(self, basis)
            testDerivatives(self, basis, deg=max(basis.getDegree() - 1, 0))
    
    def testWaveletBasis(self):
        """Test Wavelet Noboundary basis."""
        basis = pysgpp.SWaveletBase()
        testDerivatives(self, basis, max_discontinuities_count=2)
    
    def testWaveletBoundaryBasis(self):
        """Test Wavelet Boundary basis."""
        basis = pysgpp.SWaveletBoundaryBase()
        testDerivatives(self, basis, max_discontinuities_count=2)
    
    def testModWaveletBasis(self):
        """Test modified Wavelet basis."""
        basis = pysgpp.SModWaveletBase()
        testDerivatives(self, basis, max_discontinuities_count=2)
