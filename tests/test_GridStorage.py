# This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems.
#
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de)
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

class TestGridStorage(unittest.TestCase):
    def testCreateDestroy(self):
        """Test creation and destruction via storage"""
        from pysgpp import GridIndex, GridStorage
        
        i = GridIndex(1)
        s = GridStorage(1)
        
        i.set(0,1,1)
        i2 = s.create(i)
        
        self.failUnlessEqual(i.get(0), i2.get(0))
        
        s.destroy(i2)
        
    def testSerialize(self):
        """Tests serialization"""
        from pysgpp import GridStorage, HashGenerator
        
        s = GridStorage(2)
        g = HashGenerator()
        
        g.regular(s, 2)        
        
        str = s.serialize()
        
        self.failIfEqual(len(str), 0)
        
        s2 = GridStorage(str)
        
        self.failUnlessEqual(s.size(), s2.size())
        
    def testInsert(self):
        """Tests insertion of an index into a storage"""
        from pysgpp import GridIndex, GridStorage
        
        i = GridIndex(1)
        s = GridStorage(1)
        
        i.set(0,1,1)
        i2 = s.insert(i)
        
        self.failUnless(i2 == 0)
        self.failUnless(s.size() == 1)
 
    def testChilds(self):
        """Tests child construction. Superseded by GridIterator"""
        from pysgpp import GridIndex, GridStorage, HashGenerator

        s = GridStorage(1)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        i = GridIndex(1)
        
        i.set(0,1,1)
        s.left_child(i, 0)
        self.failUnlessEqual(i.get(0), [2,1])
        
        i.set(0,1,1)
        s.right_child(i, 0)
        self.failUnlessEqual(i.get(0), [2,3])
        
        
    def testTop(self):
        """Tests top construction. Now superseded by GridIterator"""
        from pysgpp import GridIndex, GridStorage, HashGenerator
        s = GridStorage(1)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        i = GridIndex(1)
        
        i.set(0,1,1)
        s.left_child(i, 0)

        s.top(i, 0)
        
        self.failUnlessEqual(i.get(0), [1,1])
        
        
    def testSeq(self):
        """Tests sequence numbers"""
        from pysgpp import GridIndex, GridStorage, HashGenerator
        s = GridStorage(1)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        i = GridIndex(1)
        
        i.set(0,1,1)
        s.left_child(i, 0)

        seq = s.seq(i)
        self.failIf(s.end(seq))

        s.left_child(i, 0)

        seq = s.seq(i)
        self.failUnless(s.end(seq))

 
class TestHashGenerator(unittest.TestCase):
    def testRegular1D(self):
        """Tests 1D grid generation"""
        from pysgpp import GridStorage, HashGenerator
        
        s = GridStorage(1)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        self.failUnlessEqual(s.size(), 3)
        
    def testRegular2D(self):
        """Tests 2D grid generation"""
        from pysgpp import GridStorage, HashGenerator
        
        s = GridStorage(2)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        self.failUnlessEqual(s.size(), 5)

        s = GridStorage(2)
        g.regular(s, 3)
        self.failUnlessEqual(s.size(), 17)
        
        from pysgpp import GridIndex
        
        i = GridIndex(2)
        i.set(0,1,1)
        i.set(1,1,1)
        
        self.failUnless(s.has_key(i))
        
        i.set(1,2,1)
        self.failUnless(s.has_key(i))
        
        i.set(0,2,3)
        self.failUnless(s.has_key(i))
        
        i.set(0,3,5)
        self.failIf(s.has_key(i))
        
        i.set(1,1,1)
        self.failUnless(s.has_key(i))

    def testRegular3D(self):
        """Tests 3D grid generation"""
        from pysgpp import GridStorage, HashGenerator
        
        s = GridStorage(3)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        self.failUnlessEqual(s.size(), 7)

class TestHashRefinement(unittest.TestCase):
    def testFreeRefine(self):
        """Tests surplus based refine for Hash-Storage"""
        from pysgpp import GridStorage, HashGenerator
        from pysgpp import SurplusRefinementFunctor, HashRefinement, DataVector
        
        s = GridStorage(2)
        g = HashGenerator()
        
        g.regular(s, 1)
        
        d = DataVector(1)
        d[0] = 1.0
        
        f = SurplusRefinementFunctor(d)
        r = HashRefinement()
        
        r.free_refine(s, f)
        
        self.failUnlessEqual(s.size(), 5)

    def testSurplusFunctor(self):
        """Tests if surplus functor correctly considers absolute values"""
        from pysgpp import GridStorage
        from pysgpp import SurplusRefinementFunctor, DataVector

        s = GridStorage(2)
        d = DataVector(1)
        f = SurplusRefinementFunctor(d)
        
        d[0] = -10.0
        self.failUnless(f(s, 0) > f.start())

        d[0] = 10.0
        self.failUnless(f(s, 0) > f.start())
        
        