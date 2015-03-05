# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


import unittest
import math

##
# Test base classes
class TestBase(unittest.TestCase):

  ##
  # General function for testing bases
  def baseTest(self, b, points):
    from pysgpp import cvar
    for t in points:
      val = b.eval(t[0], t[1], t[2])
      errorMessage = "%f != %f => (%d, %d) @ %f"%(val, t[3], t[0], t[1], t[2])
      if cvar.USING_DOUBLE_PRECISION:
        self.assertAlmostEqual(val, t[3], msg=errorMessage)
      else:
        self.assertAlmostEqual(val, t[3], msg=errorMessage, places=6)

  def linearUniformUnmodifiedTest(self, b):
    points = [(1, 1, 0.5,   1.0),
              (1, 1, 0.75,  0.5),
              (1, 1, 0.875, 0.25),
              (1, 1, 0.0,   0.0),
              (1, 1, 1.0,   0.0),
              (2, 1, 0.75,  0.0),
              (2, 3, 0.75,  1.0),
              (3, 1, 0.0,   0.0),
              (3, 1, 0.125, 1.0),
              (3, 1, 0.25,  0.0),
              ]
    self.baseTest(b, points)

  def linearModifiedTest(self, b):
    points = [(1, 1, 0.5,   1.0),
              (1, 1, 0.75,  1.0),
              (1, 1, 0.875, 1.0),
              (1, 1, 0.0,   1.0),
              (1, 1, 1.0,   1.0),
              (2, 1, 0.0,   2.0),
              (2, 1, 0.125, 1.5),
              (2, 1, 0.25,  1.0),
              (2, 1, 0.375, 0.5),
              (2, 1, 0.75,  0.0),
              (2, 3, 0.75,  1.0),
              (2, 3, 1.0,   2.0),
              (3, 1, 0.0,   2.0),
              (3, 1, 0.125, 1.0),
              (3, 1, 0.25,  0.0),
              ]
    self.baseTest(b, points)

  def linearLevelZeroTest(self, b):
    """Test boundary linear basis functions (level 0)."""
    for i in range(2):
        for x in [j/4.0 for j in range(5)]:
            self.assertEqual(b.eval(0, i, x), i*x + (1-i)*(1.0-x))

  def ccKnot(self, l, i):
    """Return Clenshaw-Curtis knot with given level and index."""
    return 0.5 * (math.cos(math.pi * (1.0 - i / 2.0**l)) + 1.0)

  def linearClenshawCurtisTest(self, b):
    points = [(1, 1, 0.5,               1.0),
              (1, 1, 0.75,              0.5),
              (1, 1, 0.875,             0.25),
              (1, 1, 0.0,               0.0),
              (1, 1, 1.0,               0.0),
              (2, 1, 0.75,              0.0),
              (2, 3, 0.75,              0.25 / (self.ccKnot(2, 3) - 0.5)),
              (2, 3, self.ccKnot(2, 3), 1.0),
              (3, 1, 0.0,               0.0),
              (3, 1, 0.125,
                   1.0 - (0.125 - self.ccKnot(3, 1)) /
                   (self.ccKnot(3, 2) - self.ccKnot(3, 1))),
              (3, 1, self.ccKnot(3, 1), 1.0),
              (3, 1, 0.25,              0.0),
              ]
    self.baseTest(b, points)

  def derivativesTest(self, b, deg=2, start_level=1, max_discontinuities_count=0):
    """Test derivatives (up to order deg, max. 2) of basis functions for level >= start_level.
    Allow for max_discontinuities_count discontinuities (e.g. for Wavelets which are cut off).
    """
    # skip test when using single precision, because then
    # the derivatives are not exact enough
    from pysgpp import cvar
    if not cvar.USING_DOUBLE_PRECISION: return
    
    dx = 1e-8
    tol = 0.01
    discontinuityTol = 1e5
    
    # levels
    for l in range(start_level, 4):
      # indices
      for i in range(1, 2**l, 2):
        f = lambda y: b.eval(l, i, y)
        
        if deg >= 1:
          # test first derivative at boundary (central difference quotient)
          df = lambda y: b.evalDx(l, i, y)
          self.assertAlmostEqual((f(2.0*dx) - f(0.0)) / (2.0*dx), df(dx), delta=tol)
          self.assertAlmostEqual((f(1.0) - f(1.0-2.0*dx)) / (2.0*dx), df(1.0-dx), delta=tol)
        if deg >= 2:
          # test second derivative at boundary (central difference quotient)
          ddf = lambda y: b.evalDxDx(l, i, y)
          self.assertAlmostEqual((df(2.0*dx) - df(0.0)) / (2.0*dx), ddf(dx), delta=tol)
          self.assertAlmostEqual((df(1.0) - df(1.0-2.0*dx)) / (2.0*dx), ddf(1.0-dx), delta=tol)
        
        # count discontinuities
        discontinuities = 0
        for x in [j/100.0 for j in range(1, 100)]:
          if abs(f(x+dx) - f(x-dx)) > discontinuityTol * dx:
            # discontinuity found
            discontinuities += 1
          else:
            # test derivatives only if the function is continuous
            if deg >= 1:
              # test first derivative (central difference quotient)
              self.assertAlmostEqual((f(x+dx) - f(x-dx)) / (2.0*dx), df(x), delta=tol)
            if deg >= 2:
              # test second derivative (central difference quotient)
              self.assertAlmostEqual((df(x+dx) - df(x-dx)) / (2.0*dx), ddf(x), delta=tol)
        self.assertLessEqual(discontinuities, max_discontinuities_count)

  def bsplinePropertiesTest(self, b, start_level=1):
    """Test basic B-spline properties (mixed monotonicity, bounds) for level >= start_level."""
    from pysgpp import cvar
    
    tol = 0.0 if cvar.USING_DOUBLE_PRECISION else 1e-4
    
    for l in range(start_level, 4):
        for i in range(1, 2**l, 2):
            # rising at the beginning
            falling = False
            fx = None
            for x in [j/100.0 for j in range(101)]:
                fx_new = b.eval(l, i, x)
                if fx is not None:
                    if falling:
                        # hope we're still falling
                        self.assertLessEqual(fx_new, fx + tol)
                    elif fx_new < fx - tol:
                        # we're now falling (and weren't until now)
                        falling = True
                fx = fx_new
                
                # test bounds
                if fx < 0.0 - tol: print b.getDegree(), l, i, x, fx
                self.assertGreaterEqual(fx, 0.0 - tol)
                
                if 1 < i < 2**l-1:
                    self.assertLessEqual(fx, 1.0 + tol)
                else:
                    # allow 2 as upper bound for leftmost and rightmost basis functions
                    # (in case of modified B-splines)
                    self.assertLessEqual(fx, 2.0 + 1e-2 + tol)

  def testLinear(self):
    from pysgpp import SLinearBase
    b = SLinearBase()
    self.linearUniformUnmodifiedTest(b)
    self.derivativesTest(b, deg=0)
  
  def testLinearBoundary(self):
    from pysgpp import SLinearBoundaryBase
    b = SLinearBoundaryBase()
    self.linearLevelZeroTest(b)
    self.linearUniformUnmodifiedTest(b)
    self.derivativesTest(b, deg=0, start_level=0)
  
  def testLinearClenshawCurtis(self):
    from pysgpp import SLinearClenshawCurtisBase
    b = SLinearClenshawCurtisBase()
    self.linearLevelZeroTest(b)
    self.linearClenshawCurtisTest(b)
    self.derivativesTest(b, deg=0, start_level=0)
  
  def testLinearModified(self):
    from pysgpp import SLinearModifiedBase
    b = SLinearModifiedBase()
    self.linearModifiedTest(b)
    self.derivativesTest(b, deg=0)

  def testLinearStretched(self):
    from pysgpp import SLinearStretchedBase
    
    b = SLinearStretchedBase()
    
# Point data format: C, L, R, value 
# (1,1) and (2,1) tried
    points = [(1.312631659, 0.5, 7, 0.125020255230769),
  (0.96716821, 0.5, 7, 0.071872032307692)
              ]
    
    self.baseTest(b, points)
        
  def testBspline(self):
    """Test B-spline Noboundary basis."""
    from pysgpp import SBsplineBase, cvar
    b = SBsplineBase(1)
    self.linearUniformUnmodifiedTest(b)
    p_max = 9 if cvar.USING_DOUBLE_PRECISION else 6
    
    for p in range(p_max+1):
      b = SBsplineBase(p)
      self.bsplinePropertiesTest(b)
      self.derivativesTest(b, deg=max(b.getDegree() - 1, 0))

  def testBsplineBoundaryBasis(self):
    """Test B-spline Boundary basis."""
    from pysgpp import SBsplineBoundaryBase, cvar
    b = SBsplineBoundaryBase(1)
    self.linearLevelZeroTest(b)
    self.linearUniformUnmodifiedTest(b)
    p_max = 9 if cvar.USING_DOUBLE_PRECISION else 6
    
    for p in range(p_max+1):
      b = SBsplineBoundaryBase(p)
      self.bsplinePropertiesTest(b, 0)
      self.derivativesTest(b, deg=max(b.getDegree() - 1, 0), start_level=0)

  def testBsplineClenshawCurtisBasis(self):
    """Test B-spline Clenshaw-Curtis basis."""
    from pysgpp import SBsplineClenshawCurtisBase
    b = SBsplineClenshawCurtisBase(1)
    self.linearLevelZeroTest(b)
    self.linearClenshawCurtisTest(b)
    
    for p in range(10):
      b = SBsplineClenshawCurtisBase(p)
      self.bsplinePropertiesTest(b, 0)
      self.derivativesTest(b, deg=max(b.getDegree() - 1, 0), start_level=0)
    
  def testBsplineModifiedBasis(self):
    """Test modified B-spline basis."""
    from pysgpp import SBsplineModifiedBase, cvar
    b = SBsplineModifiedBase(1)
    self.linearModifiedTest(b)
    p_max = 9 if cvar.USING_DOUBLE_PRECISION else 6
    
    for p in range(p_max+1):
      b = SBsplineModifiedBase(p)
      self.bsplinePropertiesTest(b)
      self.derivativesTest(b, deg=max(b.getDegree() - 1, 0))
    
  def testWaveletBasis(self):
    """Test Wavelet Noboundary basis."""
    from pysgpp import SWaveletBase
    b = SWaveletBase()
    self.derivativesTest(b, max_discontinuities_count=2)
    
  def testWaveletBoundaryBasis(self):
    """Test Wavelet Boundary basis."""
    from pysgpp import SWaveletBoundaryBase
    b = SWaveletBoundaryBase()
    self.derivativesTest(b, max_discontinuities_count=2)
    
  def testWaveletModifiedBasis(self):
    """Test modified Wavelet basis."""
    from pysgpp import SWaveletModifiedBase
    b = SWaveletModifiedBase()
    self.derivativesTest(b, max_discontinuities_count=2)

#    def testPoly(self):
#        from pysgpp import SPolyBase
#        
#        self.failUnlessRaises(Exception, SPolyBase, 0)
#        
#        b = SPolyBase(2)
#        
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 0.75),
#                  (1, 1, 0.75, 0.75),
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.125, 0.75),
#                  (2, 1, 0.25+0.125, 0.75),
#                  ]
# 
#        self.baseTest(b, points)
#        
#        b = SPolyBase(3)
#
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 0.75),
#                  (1, 1, 0.75, 0.75),
#                  
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.125, 0.875),
#                  (2, 1, 0.25+0.125, 0.625),
#                  
#                  (3, 1, 0.0625, 0.875),
#                  (3, 1, 0.125+0.0625, 0.625),
#                  
#                  (3, 3, 0.375-0.0625, 0.625),
#                  (3, 3, 0.375+0.0625, 0.875),
#                  ]
# 
#        self.baseTest(b, points)
#
#    def testModPoly(self):
#        from pysgpp import SPolyModifiedBase
#        
#        self.failUnlessRaises(Exception, SPolyModifiedBase, -1)
#        
#        b = SPolyModifiedBase(0)
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 1.0),
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.125, 1.0),
#                  (2, 3, 0.75, 1.0),
#                  (3, 1, 0.125, 1.0),
#                  ]
#
#        self.baseTest(b, points)
#
#        b = SPolyModifiedBase(1)
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 1.0),
#                  
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.0, 2.0),
#                  (2, 1, 0.5, 0.0),
#                  
#                  (2, 3, 0.75, 1.0),
#                  (2, 3, 1.0, 2.0),
#                  (2, 3, 0.5, 0.0),
#
#                  (3, 1, 0.125, 1.0),
#                  (3, 1, 0.0, 2.0),
#                  (3, 1, 0.25, 0.0),
#                  
#                  (3, 3, 0.25 + 0.125, 1.0),
#                  (3, 3, 0.5, 2.0),
#                  (3, 3, 0.25, 0.0),
#                  ]
#
#        self.baseTest(b, points)
#
#        b = SPolyModifiedBase(2)
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 1.0),
#                  
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.0, 2.0),
#                  (2, 1, 0.5, 0.0),
#                  
#                  (2, 3, 0.75, 1.0),
#                  (2, 3, 1.0, 2.0),
#                  (2, 3, 0.5, 0.0),
#                  
#                  (3, 1, 0.125, 1.0),
#                  (3, 1, 0.0, 2.0 + 2.0/3.0),
#                  (3, 1, 0.25, 0.0),
#
#                  (3, 3, 0.375, 1.0),
#                  (3, 3, 0.25, 0.0),
#                  (3, 3, 0.5, 0.0),
#                  (3, 3, (0.25+0.375)/2, 0.75)
#                  
#                  ]
#
#        self.baseTest(b, points)

        
class TestFunctions(unittest.TestCase):
  def testGetAffected(self):
    from pysgpp import HashGridIndex, HashGridStorage, SLinearBase
    from pysgpp import SGetAffectedBasisFunctions
    
    i = HashGridIndex(1)
    s = HashGridStorage(1)
    
    b = SLinearBase()
    
    i.set(0,1,1)
    s.insert(i)
    
    ga = SGetAffectedBasisFunctions(s)
    
    x = ga(b, [0.25])
    
    self.failUnlessEqual(x[0][0], 0)
    self.failUnlessEqual(x[0][1], 0.5)
        
        
  def testGetAffectedBoundaries(self):
    from pysgpp import HashGridIndex, HashGridStorage, SLinearBoundaryBase
    from pysgpp import SGetAffectedBasisFunctionsBoundaries
    
    i = HashGridIndex(1)
    s = HashGridStorage(1)
    
    b = SLinearBoundaryBase()

    i.set(0,0,0)
    s.insert(i)
    i.set(0,0,1)
    s.insert(i)        
    i.set(0,1,1)
    s.insert(i)
    
    ga = SGetAffectedBasisFunctionsBoundaries(s)
    
    x = ga(b, [0.5])
    
    self.failUnlessEqual(x[0][0], 0)
    self.failUnlessEqual(x[0][1], 0.5)
    self.failUnlessEqual(x[1][0], 1)
    self.failUnlessEqual(x[1][1], 0.5)
    self.failUnlessEqual(x[2][0], 2)
    self.failUnlessEqual(x[2][1], 1.0)   

  def testGetAffectedLinearStretched(self):
    from pysgpp import HashGridIndex, HashGridStorage, SLinearStretchedBoundaryBase
    from pysgpp import SGetAffectedBasisFunctionsLinearStretchedBoundaries
    from pysgpp import Stretching, Stretching1D, DimensionBoundary
    from pysgpp import cvar
    #print "blubb"
    str1d = Stretching1D()
    str1d.type='log'
    str1d.x_0=1
    str1d.xsi=10
    dimBound = DimensionBoundary() 
    dimBound.leftBoundary=0.5
    dimBound.rightBoundary=7
    stretch=Stretching(1,dimBound,str1d)
    i = HashGridIndex(1)
    s = HashGridStorage(1)
    s.setStretching(stretch)
      
    b = SLinearStretchedBoundaryBase()
      
    i.set(0,0,0)
    s.insert(i)
    i.set(0,0,1)
    s.insert(i)        
    i.set(0,1,1)
    s.insert(i)
    
    ga = SGetAffectedBasisFunctionsLinearStretchedBoundaries(s)
      
    x = ga(b, [0.25])
      
    self.assertEqual(x[0][0], 0)
    if cvar.USING_DOUBLE_PRECISION:
        self.assertEqual(x[0][1], 1.0384615384615385)
    else:
        self.assertAlmostEqual(x[0][1], 1.0384615384615385, places=7)

# Run tests for this file if executed as application 
if __name__=='__main__':
  unittest.main()
    
