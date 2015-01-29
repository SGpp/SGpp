# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest

class TestSHashGridIndex(unittest.TestCase):
    def testConstructor(self):
        """Tests copy constructor"""
        from pysgpp.base import HashGridIndex
        
        s = HashGridIndex(2)
        s.set(0,1,1)
        s.set(0,2,3)
        
        s2 = HashGridIndex(s)
        
        self.failUnlessEqual(s.getLevel(0), s2.getLevel(0))
        self.failUnlessEqual(s.getIndex(0), s2.getIndex(0))
        
    def testAssign(self):
        """Tests assignment"""
        from pysgpp.base import HashGridIndex
         
        s = HashGridIndex(2)
        s.set(0,2,1)
        s.set(1,2,3)
         
        s2 = HashGridIndex(5)
         
        s2.assign(s)
        self.failUnlessEqual(s.getLevel(0), s2.getLevel(0))
        self.failUnlessEqual(s.getIndex(1), s2.getIndex(1))
        
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()