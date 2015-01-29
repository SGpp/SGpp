# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org



import unittest

class TestHashGridStorage(unittest.TestCase):
    def testCreateDestroy(self):
        """Test creation and destruction via storage"""
        from pysgpp.base import HashGridIndex, HashGridStorage
        
        i = HashGridIndex(1)
        s = HashGridStorage(1)
        
        i.set(0,1,1)
        i2 = s.create(i)
        
        self.failUnlessEqual(i.get(0), i2.get(0))
        
        s.destroy(i2)
        
    def testSerialize(self):
        """Tests serialization"""
        from pysgpp.base import HashGridStorage, HashGenerator
        
        s = HashGridStorage(2)
        g = HashGenerator()
        
        g.regular(s, 2)        
        
        str = s.serialize()
        
        self.failIfEqual(len(str), 0)
        
        s2 = HashGridStorage(str)
        
        self.failUnlessEqual(s.size(), s2.size())

    def testSerializeWithLeaf(self):
        """Tests serialization with Leaf"""
        from pysgpp.base import HashGridStorage, HashGenerator
        
        srcLeaf = []
        s = HashGridStorage(2)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        for i in xrange(s.size()):
            srcLeaf.append(s.get(i).isLeaf())
        
        str = s.serialize()
        
        self.failIfEqual(len(str), 0)
        
        s2 = HashGridStorage(str)

        self.failUnlessEqual(s.size(), s2.size())
        
        for i in xrange(s.size()):
            self.failUnlessEqual(s2.get(i).isLeaf(), srcLeaf[i])        
        
        
    def testInsert(self):
        """Tests insertion of an index into a storage"""
        from pysgpp.base import HashGridIndex, HashGridStorage
        
        i = HashGridIndex(1)
        s = HashGridStorage(1)
        
        i.set(0,1,1)
        i2 = s.insert(i)
        
        self.failUnless(i2 == 0)
        self.failUnless(s.size() == 1)
 
    def testChilds(self):
        """Tests child construction. Superseded by GridIterator"""
        from pysgpp.base import HashGridIndex, HashGridStorage, HashGenerator

        s = HashGridStorage(1)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        i = HashGridIndex(1)
        
        i.set(0,1,1)
        s.left_child(i, 0)
        self.failUnlessEqual(i.get(0), [2,1])
        
        i.set(0,1,1)
        s.right_child(i, 0)
        self.failUnlessEqual(i.get(0), [2,3])
        
    def testLevelZero(self):
        """Tests child construction. Superseded by GridIterator"""
        from pysgpp.base import HashGridIndex, HashGridStorage, HashGenerator

        s = HashGridStorage(1)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        i = HashGridIndex(1)
        
        i.set(0,1,1)
        s.left_levelzero(i, 0)
        self.failUnlessEqual(i.get(0), [0,0])
        
        i.set(0,1,1)
        s.right_levelzero(i, 0)
        self.failUnlessEqual(i.get(0), [0,1])        
        
        
    def testTop(self):
        """Tests top construction. Now superseded by GridIterator"""
        from pysgpp.base import HashGridIndex, HashGridStorage, HashGenerator
        s = HashGridStorage(1)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        i = HashGridIndex(1)
        
        i.set(0,1,1)
        s.left_child(i, 0)

        s.top(i, 0)
        
        self.failUnlessEqual(i.get(0), [1,1])
        
        
    def testSeq(self):
        """Tests sequence numbers"""
        from pysgpp.base import HashGridIndex, HashGridStorage, HashGenerator
        s = HashGridStorage(1)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        i = HashGridIndex(1)
        
        i.set(0,1,1)
        s.left_child(i, 0)

        seq = s.seq(i)
        self.failIf(s.end(seq))

        s.left_child(i, 0)

        seq = s.seq(i)
        self.failUnless(s.end(seq))

 
class TestHashGenerator(unittest.TestCase):
    def testPeriodic1D(self):
        """Test 1D grid with periodic boundaries generation"""
        from pysgpp.base import HashGridStorage, HashGenerator
         
        s = HashGridStorage(1)
        g = HashGenerator()
         
        g.regularWithPeriodicBoundaries(s, 2)
         
        self.failUnlessEqual(s.size(), 4)
         
     
    def testPeriodic2D(self):
        """Tests 2D grid with periodic boundaries generation"""
        from pysgpp.base import HashGridStorage, HashGenerator
        from pysgpp.base import HashGridIndex
         
        s = HashGridStorage(2)
        g = HashGenerator()
         
        g.regularWithPeriodicBoundaries(s, 2)
        
        self.failUnlessEqual(s.size(), 12)
 
        s = HashGridStorage(2)
        g.regularWithPeriodicBoundaries(s, 3)
        self.failUnlessEqual(s.size(), 32)
        
        i = HashGridIndex(2)
        i.set(0,0,0)
        i.set(1,0,0)
        self.failUnless(s.has_key(i))
         
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
         
         
    def testPeriodic3D(self):
        """Tests 3D grid generation"""
        from pysgpp.base import HashGridStorage, HashGenerator
         
        s = HashGridStorage(3)
        g = HashGenerator()
         
        g.regularWithPeriodicBoundaries(s, 2)
         
        self.failUnlessEqual(s.size(), 32)
            
        
    def testRegular1D(self):
        """Tests 1D grid generation"""
        from pysgpp.base import HashGridStorage, HashGenerator
        
        s = HashGridStorage(1)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        self.failUnlessEqual(s.size(), 3)
        
    def testRegular2D(self):
        """Tests 2D grid generation"""
        from pysgpp.base import HashGridStorage, HashGenerator
        
        s = HashGridStorage(2)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        self.failUnlessEqual(s.size(), 5)

        s = HashGridStorage(2)
        g.regular(s, 3)
        self.failUnlessEqual(s.size(), 17)
        
        from pysgpp.base import HashGridIndex
        
        i = HashGridIndex(2)
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
        from pysgpp.base import HashGridStorage, HashGenerator
        
        s = HashGridStorage(3)
        g = HashGenerator()
        
        g.regular(s, 2)
        
        self.failUnlessEqual(s.size(), 7)
        
        
    def testRegularTruncatedBoundaries1D(self):
        """Tests 1D grid generation"""
        from pysgpp.base import HashGridStorage, HashGenerator
        
        s = HashGridStorage(1)
        g = HashGenerator()
        
        g.regularWithBoundaries(s, 2, True)
        
        self.failUnlessEqual(s.size(), 5)
        
    def testRegularTruncatedBoundaries2D(self):
        """Tests 2D grid generation"""
        from pysgpp.base import HashGridStorage, HashGenerator
        
        s = HashGridStorage(2)
        g = HashGenerator()
        
        g.regularWithBoundaries(s, 2, True)
        
        self.failUnlessEqual(s.size(), 21)

        s = HashGridStorage(2)
        g.regularWithBoundaries(s, 3, True)
        self.failUnlessEqual(s.size(), 49)
        
        from pysgpp.base import HashGridIndex
        
        i = HashGridIndex(2)
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
        
        i.set(1,0,0)
        self.failUnless(s.has_key(i))
        
        i.set(0,0,0)
        self.failUnless(s.has_key(i))


    def testRegularTruncatedBoundaries3D(self):
        """Tests 3D grid generation"""
        from pysgpp.base import HashGridStorage, HashGenerator
        
        s = HashGridStorage(3)
        g = HashGenerator()
        
        g.regularWithBoundaries(s, 2, True)
        
        self.failUnlessEqual(s.size(), 81)        



class TestHashRefinement(unittest.TestCase):
    def testFreeRefine(self):
        """Tests surplus based refine for Hash-Storage"""
        from pysgpp.base import HashGridStorage, HashGenerator
        from pysgpp.base import SurplusRefinementFunctor, HashRefinement, DataVector
        
        s = HashGridStorage(2)
        g = HashGenerator()
        
        g.regular(s, 1)
        
        d = DataVector(1)
        d[0] = 1.0
        
        f = SurplusRefinementFunctor(d)
        r = HashRefinement()
        
        r.free_refine(s, f)
        
        self.failUnlessEqual(s.size(), 5)
        
    def testFreeRefineTruncatedBoundaries(self):
        """Tests surplus based refine for Hash-Storage"""
        from pysgpp.base import HashGridStorage, HashGenerator
        from pysgpp.base import SurplusRefinementFunctor, HashRefinementBoundaries, DataVector
        
        s = HashGridStorage(2)
        g = HashGenerator()
        
        g.regularWithBoundaries(s, 1, True)
        
        d = DataVector(9)
        d[0] = 0.0
        d[1] = 0.0
        d[2] = 0.0
        d[3] = 0.0
        d[4] = 0.0
        d[5] = 0.0
        d[6] = 0.0
        d[7] = 0.0
        d[8] = 1.0
        
        f = SurplusRefinementFunctor(d)
        r = HashRefinementBoundaries()
        
        r.free_refine(s, f)
        
        self.failUnlessEqual(s.size(), 21)     
        
    def testFreeRefineTruncatedBoundaries(self):
        """Tests surplus based refine for Hash-Storage"""
        from pysgpp.base import HashGridStorage, HashGenerator
        from pysgpp.base import SurplusRefinementFunctor, HashRefinementBoundaries, DataVector
        
        s = HashGridStorage(2)
        g = HashGenerator()
        
        g.regularWithBoundaries(s, 2, False)
        
        d = DataVector(17)

        for i in xrange(len(d)):
            d[i] = 0.0
            
        d[12] = 1.0
        
        f = SurplusRefinementFunctor(d)
        r = HashRefinementBoundaries()
        
        r.free_refine(s, f)
        
        self.failUnlessEqual(s.size(), 21)          
           

    def testSurplusFunctor(self):
        """Tests if surplus functor correctly considers absolute values"""
        from pysgpp.base import HashGridStorage
        from pysgpp.base import SurplusRefinementFunctor, DataVector

        s = HashGridStorage(2)
        d = DataVector(1)
        f = SurplusRefinementFunctor(d)
        
        d[0] = -10.0
        self.failUnless(f(s, 0) > f.start())

        d[0] = 10.0
        self.failUnless(f(s, 0) > f.start())
        
if __name__=="__main__":
    unittest.main()         
