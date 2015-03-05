# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
from pysgpp import *

class TestGridFactory(unittest.TestCase):
    def testCreation(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(2)
        self.failIfEqual(factory, None)
        
        storage = factory.getStorage()
        self.failIfEqual(storage, None)
        
    def testSerializationLinear(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
    def testSerializationModLinear(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createModLinearGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
    def testSerializationLinearTruncatedBoundary(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createLinearTruncatedBoundaryGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
    def testSerializationLinearBoundary(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
    def testSerializationPrewavelet(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createPrewaveletGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
#    def testSerializationPoly(self):
#        from pysgpp import Grid
#        
#        factory = Grid.createPolyGrid(2,2)
#        self.failIfEqual(factory, None)
#
#        gen = factory.createGridGenerator()
#        gen.regular(3)
#
#        str = factory.serialize()
#        self.assert_(len(str) > 0)
#        
#        newfac = Grid.unserialize(str)
#        self.failIfEqual(newfac, None)
#        
#        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())

    def testSerializationLinearBoundingBox(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)
        
        boundingBox = factory.getBoundingBox()
        tempBound = boundingBox.getBoundary(0)
        tempBound.leftBoundary = 0.0
        tempBound.rightBoundary = 100.0
        tempBound.bDirichletLeft = False;
        tempBound.bDirichletRight = False;
        boundingBox.setBoundary(0, tempBound)
        
        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
        boundingBox = newfac.getBoundingBox()    
        tempBound = boundingBox.getBoundary(0)
        self.assertEqual(0.0, tempBound.leftBoundary)
        self.assertEqual(100.0, tempBound.rightBoundary)
        self.assertEqual(False, tempBound.bDirichletLeft)
        self.assertEqual(False, tempBound.bDirichletRight)
        
    def testSerializationModLinearBoundingBox(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createModLinearGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        boundingBox = factory.getBoundingBox()
        tempBound = boundingBox.getBoundary(0)
        tempBound.leftBoundary = 0.0
        tempBound.rightBoundary = 100.0
        tempBound.bDirichletLeft = False;
        tempBound.bDirichletRight = False;
        boundingBox.setBoundary(0, tempBound)

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
        boundingBox = newfac.getBoundingBox()
        tempBound = boundingBox.getBoundary(0)
        self.assertEqual(0.0, tempBound.leftBoundary)
        self.assertEqual(100.0, tempBound.rightBoundary)
        self.assertEqual(False, tempBound.bDirichletLeft)
        self.assertEqual(False, tempBound.bDirichletRight)
        
    def testSerializationLinearTruncatedBoundaryBoundingBox(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createLinearTruncatedBoundaryGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        boundingBox = factory.getBoundingBox()
        tempBound = boundingBox.getBoundary(0)
        tempBound.leftBoundary = 0.0
        tempBound.rightBoundary = 100.0
        tempBound.bDirichletLeft = False;
        tempBound.bDirichletRight = False;
        boundingBox.setBoundary(0, tempBound)

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
        boundingBox = newfac.getBoundingBox()
        tempBound = boundingBox.getBoundary(0)
        self.assertEqual(0.0, tempBound.leftBoundary)
        self.assertEqual(100.0, tempBound.rightBoundary)
        self.assertEqual(False, tempBound.bDirichletLeft)
        self.assertEqual(False, tempBound.bDirichletRight)
        
    def testSerializationLinearBoundaryBoundingBox(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        boundingBox = factory.getBoundingBox()
        tempBound = boundingBox.getBoundary(0)
        tempBound.leftBoundary = 0.0
        tempBound.rightBoundary = 100.0
        tempBound.bDirichletLeft = False;
        tempBound.bDirichletRight = False;
        boundingBox.setBoundary(0, tempBound)

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
            
        boundingBox = newfac.getBoundingBox()
        tempBound = boundingBox.getBoundary(0)
        self.assertEqual(0.0, tempBound.leftBoundary)
        self.assertEqual(100.0, tempBound.rightBoundary)
        self.assertEqual(False, tempBound.bDirichletLeft)
        self.assertEqual(False, tempBound.bDirichletRight)
        
#    def testSerializationPolyBoundingBox(self):
#        from pysgpp import Grid
#        
#        factory = Grid.createPolyGrid(2,2)
#        self.failIfEqual(factory, None)
#
#        gen = factory.createGridGenerator()
#        gen.regular(3)
#
#        boundingBox = factory.getBoundingBox()
#        tempBound = boundingBox.getBoundary(0)
#        tempBound.leftBoundary = 0.0
#        tempBound.rightBoundary = 100.0
#        tempBound.bDirichletLeft = False;
#        tempBound.bDirichletRight = False;
#        boundingBox.setBoundary(0, tempBound)
#
#        str = factory.serialize()
#        self.assert_(len(str) > 0)
#        
#        newfac = Grid.unserialize(str)
#        self.failIfEqual(newfac, None)
#        
#        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
#        
#        boundingBox = newfac.getBoundingBox()
#        tempBound = boundingBox.getBoundary(0)
#        self.assertEqual(0.0, tempBound.leftBoundary)
#        self.assertEqual(100.0, tempBound.rightBoundary)
#        self.assertEqual(False, tempBound.bDirichletLeft)
#        self.assertEqual(False, tempBound.bDirichletRight)

    def testSerializationLinearWithLeaf(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        srcLeaf = []
        factory = Grid.createLinearGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        for i in xrange(factory.getStorage().size()):
            srcLeaf.append(factory.getStorage().get(i).isLeaf())

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
        for i in xrange(factory.getStorage().size()):
            self.failUnlessEqual(newfac.getStorage().get(i).isLeaf(), srcLeaf[i])
        
    def testSerializationModLinearWithLeaf(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        srcLeaf = []
        factory = Grid.createModLinearGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        for i in xrange(factory.getStorage().size()):
            srcLeaf.append(factory.getStorage().get(i).isLeaf())

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
        for i in xrange(factory.getStorage().size()):
            self.failUnlessEqual(newfac.getStorage().get(i).isLeaf(), srcLeaf[i])
        
    def testSerializationLinearTruncatedBoundaryWithLeaf(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        srcLeaf = []
        factory = Grid.createLinearTruncatedBoundaryGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        for i in xrange(factory.getStorage().size()):
            srcLeaf.append(factory.getStorage().get(i).isLeaf())

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
        for i in xrange(factory.getStorage().size()):
            self.failUnlessEqual(newfac.getStorage().get(i).isLeaf(), srcLeaf[i])
        
    def testSerializationLinearBoundaryWithLeaf(self):
        """Uses Linear grid for tests"""
        from pysgpp import Grid
        
        srcLeaf = []
        factory = Grid.createLinearBoundaryGrid(2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        for i in xrange(factory.getStorage().size()):
            srcLeaf.append(factory.getStorage().get(i).isLeaf())

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
        for i in xrange(factory.getStorage().size()):
            self.failUnlessEqual(newfac.getStorage().get(i).isLeaf(), srcLeaf[i])       
        
#    def testSerializationPolyWithLeaf(self):
#        from pysgpp import Grid
#        
#        srcLeaf = []
#        factory = Grid.createPolyGrid(2,2)
#        self.failIfEqual(factory, None)
#
#        gen = factory.createGridGenerator()
#        gen.regular(3)
#        
#        for i in xrange(factory.getStorage().size()):
#            srcLeaf.append(factory.getStorage().get(i).isLeaf())
#
#        str = factory.serialize()
#        self.assert_(len(str) > 0)
#        
#        newfac = Grid.unserialize(str)
#        self.failIfEqual(newfac, None)
#        
#        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
#        
#        for i in xrange(factory.getStorage().size()):
#            self.failUnlessEqual(newfac.getStorage().get(i).isLeaf(), srcLeaf[i])
        
#class TestPolyGrid(unittest.TestCase):
#    def testGeneration(self):
#        from pysgpp import Grid
#
#        factory = Grid.createPolyGrid(2, 4)
#        self.failIfEqual(factory, None)
#        
#        storage = factory.getStorage()
#        self.failIfEqual(storage, None)
#        
#        op = factory.createOperationB()
#        self.failIfEqual(op, None)
#
#        op = factory.createOperationEval()
#        self.failIfEqual(op, None)
#        
#        self.failUnlessRaises(Exception, factory.createOperationLaplace)
#
#
#class TestModPolyGrid(unittest.TestCase):
#    def testGeneration(self):
#        from pysgpp import Grid
#
#        factory = Grid.createModPolyGrid(2, 4)
#        self.failIfEqual(factory, None)
#        
#        storage = factory.getStorage()
#        self.failIfEqual(storage, None)
#        
#        op = factory.createOperationB()
#        self.failIfEqual(op, None)
#
#        op = factory.createOperationEval()
#        self.failIfEqual(op, None)
#        
#        self.failUnlessRaises(Exception, factory.createOperationLaplace)
#
#
#class TestModLinearGrid(unittest.TestCase):
#    def testGeneration(self):
#        from pysgpp import Grid
#
#        factory = Grid.createModLinearGrid(2)
#        self.failIfEqual(factory, None)
#        
#        storage = factory.getStorage()
#        self.failIfEqual(storage, None)
#        
#        op = factory.createOperationB()
#        self.failIfEqual(op, None)
#
#        op = factory.createOperationEval()
#        self.failIfEqual(op, None)
#        
#        op = factory.createOperationLaplace()
#        self.failIfEqual(op, None)
        
        
class TestLinearGrid(unittest.TestCase):
    def testGeneration(self):
        from pysgpp import Grid, DataVector
        factory = Grid.createLinearGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        self.failIfEqual(gen, None)
        
        self.failUnlessEqual(storage.size(), 0)
        gen.regular(3)
        self.failUnlessEqual(storage.size(), 17)
        
        #This should fail
        self.failUnlessRaises(Exception, gen.regular, 3)
        
    def testRefinement(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor
        factory = Grid.createLinearGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        self.failUnlessEqual(storage.size(), 1)
        alpha = DataVector(1)
        alpha[0] = 1.0
        func = SurplusRefinementFunctor(alpha)
        
        gen.refine(func)
        self.failUnlessEqual(storage.size(), 5)

    def testOperationMultipleEval(self):
        from pysgpp import Grid, DataVector, DataMatrix
        factory = Grid.createLinearGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(2)
        
        alpha = DataVector(factory.getStorage().size())
        p = DataMatrix(1,1)
        beta = DataVector(1)
        
        
        alpha.setAll(0.0)
        p.set(0,0,0.25)
        beta[0] = 1.0
        
        opb = createOperationMultipleEval(factory, p)
        opb.multTranspose(beta, alpha)
        
        self.failUnlessAlmostEqual(alpha[0], 0.5)
        self.failUnlessAlmostEqual(alpha[1], 1.0)
        self.failUnlessAlmostEqual(alpha[2], 0.0)
        
        alpha.setAll(0.0)
        alpha[0] = 1.0
        
        p.set(0,0, 0.25)
        
        beta[0] = 0.0
        
        opb.mult(alpha, beta)
        self.failUnlessAlmostEqual(beta[0], 0.5)

    def testOperationTest_test(self):
        from pysgpp import Grid, DataVector, DataMatrix

        factory = Grid.createLinearGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        
        data = DataMatrix(1,1)
        data.setAll(0.25)
        classes = DataVector(1)
        classes.setAll(1.0)

        testOP = createOperationTest(factory)

        alpha.setAll(1.0)
        c = testOP.test(alpha, data, classes)
        self.failUnless(c > 0.0)
        
        alpha.setAll(-1.0)
        c = testOP.test(alpha, data, classes)
        self.failUnless(c == 0.0)
        
    def testOperationEval_eval(self):
        from pysgpp import Grid, DataVector

        factory = Grid.createLinearGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        alpha.setAll(1.0)

        p = DataVector(1)
        p.setAll(0.25)
        
        eval = createOperationEval(factory)

        self.failUnlessAlmostEqual(eval.eval(alpha, p), 0.5)

class TestLinearTruncatedBoundaryGrid(unittest.TestCase):
    def testGeneration(self):
        from pysgpp import Grid, DataVector
        factory = Grid.createLinearTruncatedBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        self.failIfEqual(gen, None)
        
        self.failUnlessEqual(storage.size(), 0)
        gen.regular(3)
        self.failUnlessEqual(storage.size(), 49)
        
        #This should fail
        self.failUnlessRaises(Exception, gen.regular, 3)
        
    def testRefinement2d(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor
        factory = Grid.createLinearTruncatedBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        self.failUnlessEqual(storage.size(), 9)
        alpha = DataVector(9)
        alpha[0] = 0.0
        alpha[1] = 0.0
        alpha[2] = 0.0
        alpha[3] = 0.0
        alpha[4] = 0.0
        alpha[5] = 0.0
        alpha[6] = 0.0
        alpha[7] = 0.0
        alpha[8] = 1.0
        func = SurplusRefinementFunctor(alpha)
        
        gen.refine(func)
        self.failUnlessEqual(storage.size(), 21)
        
    def testRefinement3d(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor
        factory = Grid.createLinearTruncatedBoundaryGrid(3)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        self.failUnlessEqual(storage.size(), 27)
        alpha = DataVector(27)
        for i in xrange(len(alpha)):
             alpha[i] = 0.0

        alpha[26] = 1.0
        func = SurplusRefinementFunctor(alpha)
        
        gen.refine(func)
        self.failUnlessEqual(storage.size(), 81)

    def testOperationMultipleEval(self):
        from pysgpp import Grid, DataVector, DataMatrix
        factory = Grid.createLinearTruncatedBoundaryGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(2)
        
        alpha = DataVector(factory.getStorage().size())
        p = DataMatrix(1,1)
        beta = DataVector(1)
        
        
        alpha.setAll(0.0)
        p.set(0,0,0.25)
        beta[0] = 1.0
        
        opb = createOperationMultipleEval(factory, p)
        opb.multTranspose(beta, alpha)
        
        self.failUnlessAlmostEqual(alpha[0], 0.75)
        self.failUnlessAlmostEqual(alpha[1], 0.25)
        self.failUnlessAlmostEqual(alpha[2], 0.5)
        self.failUnlessAlmostEqual(alpha[3], 1.0)
        self.failUnlessAlmostEqual(alpha[4], 0.0)
        
        alpha.setAll(0.0)
        alpha[2] = 1.0
        
        p.set(0,0, 0.25)
        
        beta[0] = 0.0
        
        opb.mult(alpha, beta)
        self.failUnlessAlmostEqual(beta[0], 0.5)

    def testOperationTest_test(self):
        from pysgpp import Grid, DataVector, DataMatrix

        factory = Grid.createLinearTruncatedBoundaryGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        
        data = DataMatrix(1,1)
        data.setAll(0.25)
        classes = DataVector(1)
        classes.setAll(1.0)

        testOP = createOperationTest(factory)

        alpha[0] = 0.0
        alpha[1] = 0.0
        alpha[2] = 1.0
        
        c = testOP.test(alpha, data, classes)
        self.failUnless(c > 0.0)
        
        alpha[0] = 0.0
        alpha[1] = 0.0
        alpha[2] = -1.0
        c = testOP.test(alpha, data, classes)
        self.failUnless(c == 0.0)
        
    def testOperationEval_eval(self):
        from pysgpp import Grid, DataVector

        factory = Grid.createLinearTruncatedBoundaryGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        alpha.setAll(1.0)

        p = DataVector(1)
        p.setAll(0.25)
        
        eval = createOperationEval(factory)

        self.failUnlessAlmostEqual(eval.eval(alpha, p), 1.5)

class TestLinearBoundaryGrid(unittest.TestCase):
    def testGeneration(self):
        from pysgpp import Grid, DataVector
        factory = Grid.createLinearBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        self.failIfEqual(gen, None)
        
        self.failUnlessEqual(storage.size(), 0)
        gen.regular(3)
        self.failUnlessEqual(storage.size(), 37)
        
        #This should fail
        self.failUnlessRaises(Exception, gen.regular, 3)
        
    def testRefinement2d_one(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor
        factory = Grid.createLinearBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(0)
        
        self.failUnlessEqual(storage.size(), 4)
        
        alpha = DataVector(4)
    
        for i in xrange(len(alpha)):
            alpha[i] = 0.0

        alpha[0] = 1.0
        func = SurplusRefinementFunctor(alpha)
            
        gen.refine(func)
        self.failUnlessEqual(storage.size(), 8)
        
    def testRefinement2d_two(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor
        factory = Grid.createLinearBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(0)
        
        alpha = DataVector(4)
    
        for i in xrange(len(alpha)):
            alpha[i] = 0.0

        alpha[0] = 1.0
        func = SurplusRefinementFunctor(alpha)
            
        gen.refine(func)
        
        alpha2 = DataVector(8)
    
        for i in xrange(len(alpha2)):
            alpha2[i] = 0.0

        alpha2[4] = 1.0
        func = SurplusRefinementFunctor(alpha2)
            
        gen.refine(func)    
        self.failUnlessEqual(storage.size(), 13)        

    def testRefinement2d_three(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor
        factory = Grid.createLinearBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(0)
        
        alpha = DataVector(4)
    
        for i in xrange(len(alpha)):
            alpha[i] = 0.0

        alpha[0] = 1.0
        func = SurplusRefinementFunctor(alpha)
            
        gen.refine(func)
        
        alpha2 = DataVector(8)
    
        for i in xrange(len(alpha2)):
            alpha2[i] = 0.0

        alpha2[4] = 1.0
        func = SurplusRefinementFunctor(alpha2)
            
        gen.refine(func)
        
        alpha3 = DataVector(13)
    
        for i in xrange(len(alpha3)):
             alpha3[i] = 0.0

        alpha3[11] = 1.0
        func = SurplusRefinementFunctor(alpha3)
            
        gen.refine(func)   
        self.failUnlessEqual(storage.size(), 18) 

    def testRefinement2d_four(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor
        factory = Grid.createLinearBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(0)
        
        alpha = DataVector(4)
    
        for i in xrange(len(alpha)):
            alpha[i] = 0.0

        alpha[0] = 1.0
        func = SurplusRefinementFunctor(alpha)
            
        gen.refine(func)
        
        alpha2 = DataVector(8)
    
        for i in xrange(len(alpha2)):
            alpha2[i] = 0.0

        alpha2[4] = 1.0
        func = SurplusRefinementFunctor(alpha2)
            
        gen.refine(func)
        
        alpha3 = DataVector(13)
    
        for i in xrange(len(alpha3)):
             alpha3[i] = 0.0

        alpha3[11] = 1.0
        func = SurplusRefinementFunctor(alpha3)
            
        gen.refine(func)
        
        alpha4 = DataVector(18)
    
        for i in xrange(len(alpha4)):
            alpha4[i] = 0.0

        alpha4[12] = 1.0
        func = SurplusRefinementFunctor(alpha4)
            
        gen.refine(func)
        self.failUnlessEqual(storage.size(), 25) 

    def testRefinement2d_five(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor
        factory = Grid.createLinearBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(0)
        
        alpha = DataVector(4)
    
        for i in xrange(len(alpha)):
            alpha[i] = 0.0

        alpha[0] = 1.0
        func = SurplusRefinementFunctor(alpha)
            
        gen.refine(func)
        
        alpha2 = DataVector(8)
    
        for i in xrange(len(alpha2)):
            alpha2[i] = 0.0

        alpha2[4] = 1.0
        func = SurplusRefinementFunctor(alpha2)
            
        gen.refine(func)
        
        alpha3 = DataVector(13)
    
        for i in xrange(len(alpha3)):
             alpha3[i] = 0.0

        alpha3[11] = 1.0
        func = SurplusRefinementFunctor(alpha3)
            
        gen.refine(func)
        
        alpha4 = DataVector(18)
    
        for i in xrange(len(alpha4)):
            alpha4[i] = 0.0

        alpha4[12] = 1.0
        func = SurplusRefinementFunctor(alpha4)
            
        gen.refine(func)
        
        alpha5 = DataVector(25)
    
        for i in xrange(len(alpha5)):
            alpha5[i] = 0.0

        alpha5[23] = 1.0
        func = SurplusRefinementFunctor(alpha5)
            
        gen.refine(func)
        self.failUnlessEqual(storage.size(), 29) 

    def testOperationMultipleEval(self):
        from pysgpp import Grid, DataVector, DataMatrix
        factory = Grid.createLinearBoundaryGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(2)
        
        alpha = DataVector(factory.getStorage().size())
        p = DataMatrix(1,1)
        beta = DataVector(1)
        
        
        alpha.setAll(0.0)
        p.set(0,0,0.25)
        beta[0] = 1.0
        
        opb = createOperationMultipleEval(factory, p)
        opb.multTranspose(beta, alpha)
        
        self.failUnlessAlmostEqual(alpha[0], 0.75)
        self.failUnlessAlmostEqual(alpha[1], 0.25)
        self.failUnlessAlmostEqual(alpha[2], 0.5)
        self.failUnlessAlmostEqual(alpha[3], 1.0)
        self.failUnlessAlmostEqual(alpha[4], 0.0)
        
        alpha.setAll(0.0)
        alpha[2] = 1.0
        
        p.set(0,0, 0.25)
        
        beta[0] = 0.0
        
        opb.mult(alpha, beta)
        self.failUnlessAlmostEqual(beta[0], 0.5)

    def testOperationTest_test(self):
        from pysgpp import Grid, DataVector, DataMatrix

        factory = Grid.createLinearBoundaryGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        
        data = DataMatrix(1,1)
        data.setAll(0.25)
        classes = DataVector(1)
        classes.setAll(1.0)

        testOP = createOperationTest(factory)

        alpha[0] = 0.0
        alpha[1] = 0.0
        alpha[2] = 1.0
        
        c = testOP.test(alpha, data, classes)
        self.failUnless(c > 0.0)
        
        alpha[0] = 0.0
        alpha[1] = 0.0
        alpha[2] = -1.0
        c = testOP.test(alpha, data, classes)
        self.failUnless(c == 0.0)
        
    def testOperationEval_eval(self):
        from pysgpp import Grid, DataVector

        factory = Grid.createLinearBoundaryGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        alpha.setAll(1.0)

        p = DataVector(1)
        p.setAll(0.25)
        
        eval = createOperationEval(factory)

        self.failUnlessAlmostEqual(eval.eval(alpha, p), 1.5)

class TestLinearStretchedTruncatedBoundaryGrid(unittest.TestCase):
    def testGeneration(self):
        from pysgpp import Grid, DataVector, Stretching, Stretching1D, DimensionBoundary, DimensionBoundaryVector, Stretching1DVector 
	str1d = Stretching1D()
	str1d.type='log'
	str1d.x_0=1
	str1d.xsi=10
	dimBound = DimensionBoundary() 
 	dimBound.leftBoundary=0.5
	dimBound.rightBoundary=7
	dimBoundaryVector=DimensionBoundaryVector(2)
	dimBoundaryVector[0]=dimBound;
	dimBoundaryVector[1]=dimBound;
	str1dvector = Stretching1DVector(2)
	str1dvector[0]=str1d
	str1dvector[1]=str1d
	stretch = Stretching(2, dimBoundaryVector, str1dvector)

        factory = Grid.createLinearTruncatedBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        self.failIfEqual(gen, None)
        
        self.failUnlessEqual(storage.size(), 0)
        gen.regular(3)
        self.failUnlessEqual(storage.size(), 49)
        
        #This should fail
        self.failUnlessRaises(Exception, gen.regular, 3)
        
    def testRefinement2d(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor 
        factory = Grid.createLinearStretchedTruncatedBoundaryGrid(2)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        self.failUnlessEqual(storage.size(), 9)
        alpha = DataVector(9)
        alpha[0] = 0.0
        alpha[1] = 0.0
        alpha[2] = 0.0
        alpha[3] = 0.0
        alpha[4] = 0.0
        alpha[5] = 0.0
        alpha[6] = 0.0
        alpha[7] = 0.0
        alpha[8] = 1.0
        func = SurplusRefinementFunctor(alpha)
        
        gen.refine(func)
        self.failUnlessEqual(storage.size(), 21)
        
    def testRefinement3d(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor
        factory = Grid.createLinearStretchedTruncatedBoundaryGrid(3)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        self.failUnlessEqual(storage.size(), 27)
        alpha = DataVector(27)
        for i in xrange(len(alpha)):
             alpha[i] = 0.0

        alpha[26] = 1.0
        func = SurplusRefinementFunctor(alpha)
        
        gen.refine(func)
        self.failUnlessEqual(storage.size(), 81)

    def testOperationMultipleEval(self):
        from pysgpp import Grid, DataVector, DataMatrix, Stretching, Stretching1D, DimensionBoundary

#Create Stretching
	str1d = Stretching1D()
	str1d.type='log'
	str1d.x_0=1
	str1d.xsi=10
	dimBound = DimensionBoundary() 
 	dimBound.leftBoundary=0.5
	dimBound.rightBoundary=7
	stretch=Stretching(1,dimBound,str1d)

        factory = Grid.createLinearStretchedTruncatedBoundaryGrid(1)
	factory.getStorage().setStretching(stretch)
        gen = factory.createGridGenerator()
        gen.regular(2)
        
        alpha = DataVector(factory.getStorage().size())
        p = DataMatrix(1,1)
        beta = DataVector(1)
        
        
        alpha.setAll(0.0)
        p.set(0,0,0.25)
        beta[0] = 1.0
        
        opb = createOperationMultipleEval(factory, p)
        opb.multTranspose(beta, alpha)
        
        self.failUnlessAlmostEqual(alpha[0], 1.038461538)
        self.failUnlessAlmostEqual(alpha[1], -0.0384615)
        self.failUnlessAlmostEqual(alpha[2], -0.1823714)
        self.failUnlessAlmostEqual(alpha[3], -0.53513915)
        self.failUnlessAlmostEqual(alpha[4], 0.0)
        
        alpha.setAll(0.0)
        alpha[2] = 1.0
        
        p.set(0,0, 0.25)
        
        beta[0] = 0.0
        
        opb.mult(alpha, beta)
        self.failUnlessAlmostEqual(beta[0],-0.182371437)

    def testOperationTest_test(self):
        from pysgpp import Grid, DataVector, DataMatrix, Stretching, Stretching1D, DimensionBoundary

     	str1d = Stretching1D()
	str1d.type='log'
	str1d.x_0=1
	str1d.xsi=10
	dimBound = DimensionBoundary() 
 	dimBound.leftBoundary=0.5
	dimBound.rightBoundary=7
	stretch=Stretching(1,dimBound,str1d)

        factory = Grid.createLinearStretchedTruncatedBoundaryGrid(1)
	factory.getStorage().setStretching(stretch)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        
        data = DataMatrix(1,1)
        data.setAll(0.25)
        classes = DataVector(1)
        classes.setAll(1.0)

        testOP = createOperationTest(factory)

        alpha[0] = 0.0
        alpha[1] = 0.0
        alpha[2] = 1.0
        
        c = testOP.test(alpha, data, classes)
        #print c
        self.failUnless(c > 0.0)
        
        alpha[0] = 0.0
        alpha[1] = 0.0
        alpha[2] = -1.0
        c = testOP.test(alpha, data, classes)
        self.failUnless(c == 0.0)
        
    def testOperationEval_eval(self):
        from pysgpp import Grid, DataVector, Stretching, Stretching1D, DimensionBoundary

     	str1d = Stretching1D()
	str1d.type='log'
	str1d.x_0=1
	str1d.xsi=10
	dimBound = DimensionBoundary() 
 	dimBound.leftBoundary=0.5
	dimBound.rightBoundary=7
	stretch=Stretching(1,dimBound,str1d)

        factory = Grid.createLinearStretchedTruncatedBoundaryGrid(1)
	factory.getStorage().setStretching(stretch)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        alpha.setAll(1.0)

        p = DataVector(1)
        p.setAll(0.25)
        
        eval = createOperationEval(factory)

        self.failUnlessAlmostEqual(eval.eval(alpha, p), 0.8176285620)

#unittest.main()
