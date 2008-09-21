# This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems.
#
# Copyright (C) 2008  Joerg Blank (blankj@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with pyclass. If not, see <http://www.gnu.org/licenses/>.
#


import unittest

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
        
    def testSerializationPoly(self):
        from pysgpp import Grid
        
        factory = Grid.createPolyGrid(2,2)
        self.failIfEqual(factory, None)

        gen = factory.createGridGenerator()
        gen.regular(3)

        str = factory.serialize()
        self.assert_(len(str) > 0)
        
        newfac = Grid.unserialize(str)
        self.failIfEqual(newfac, None)
        
        self.assertEqual(factory.getStorage().size(), newfac.getStorage().size())
        
class TestPolyGrid(unittest.TestCase):
    def testGeneration(self):
        from pysgpp import Grid

        factory = Grid.createPolyGrid(2, 4)
        self.failIfEqual(factory, None)
        
        storage = factory.getStorage()
        self.failIfEqual(storage, None)
        
        op = factory.createOperationB()
        self.failIfEqual(op, None)

        op = factory.createOperationEval()
        self.failIfEqual(op, None)
        
        self.failUnlessRaises(Exception, factory.createOperationLaplace)


class TestModPolyGrid(unittest.TestCase):
    def testGeneration(self):
        from pysgpp import Grid

        factory = Grid.createModPolyGrid(2, 4)
        self.failIfEqual(factory, None)
        
        storage = factory.getStorage()
        self.failIfEqual(storage, None)
        
        op = factory.createOperationB()
        self.failIfEqual(op, None)

        op = factory.createOperationEval()
        self.failIfEqual(op, None)
        
        self.failUnlessRaises(Exception, factory.createOperationLaplace)


class TestModLinearGrid(unittest.TestCase):
    def testGeneration(self):
        from pysgpp import Grid

        factory = Grid.createModLinearGrid(2)
        self.failIfEqual(factory, None)
        
        storage = factory.getStorage()
        self.failIfEqual(storage, None)
        
        op = factory.createOperationB()
        self.failIfEqual(op, None)

        op = factory.createOperationEval()
        self.failIfEqual(op, None)
        
        op = factory.createOperationLaplace()
        self.failIfEqual(op, None)
        
        
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

    def testOperationB(self):
        from pysgpp import Grid, DataVector
        factory = Grid.createLinearGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(2)
        
        alpha = DataVector(factory.getStorage().size())
        p = DataVector(1,1)
        beta = DataVector(1)
        
        
        alpha.setAll(0.0)
        p[0] = 0.25
        beta[0] = 1.0
        
        opb = factory.createOperationB()
        opb.mult(beta, p, alpha)
        
        self.failUnlessAlmostEqual(alpha[0], 0.5)
        self.failUnlessAlmostEqual(alpha[1], 1.0)
        self.failUnlessAlmostEqual(alpha[2], 0.0)
        
        alpha.setAll(0.0)
        alpha[0] = 1.0
        
        p[0] = 0.25
        
        beta[0] = 0.0
        
        opb.multTranspose(alpha, p, beta)
        self.failUnlessAlmostEqual(beta[0], 0.5)

    def testOperationEval_test(self):
        from pysgpp import Grid, DataVector

        factory = Grid.createLinearGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        
        data = DataVector(1,1)
        data.setAll(0.25)
        classes = DataVector(1,1)
        classes.setAll(1.0)

        eval = factory.createOperationEval()

        alpha.setAll(1.0)
        c = eval.test(alpha, data, classes)
        self.failUnless(c > 0.0)
        
        alpha.setAll(-1.0)
        c = eval.test(alpha, data, classes)
        self.failUnless(c == 0.0)
        
    def testOperationEval_eval(self):
        from pysgpp import Grid, DataVector

        factory = Grid.createLinearGrid(1)
        gen = factory.createGridGenerator()
        gen.regular(1)
        
        alpha = DataVector(factory.getStorage().size())        
        alpha.setAll(1.0)

        p = DataVector(1,1)
        p.setAll(0.25)
        
        eval = factory.createOperationEval()

        self.failUnlessAlmostEqual(eval.eval(alpha, p), 0.5)


