# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest

class TestSGridIndex(unittest.TestCase):
    def testConstructor(self):
        """Tests copy constructor"""
        from pysgpp.base import GridIndex
        
        s = GridIndex(2)
        s.set(0,1,1)
        s.set(0,2,3)
        
        s2 = GridIndex(s)
        self.failUnlessEqual(s.get(0), s2.get(0))
        
    def testAssign(self):
        """Tests assignment"""
        from pysgpp.base import GridIndex
        
        s = GridIndex(2)
        s.set(0,1,1)
        s.set(1,2,3)
        
        s2 = GridIndex(5)
        
        s2.assign(s)
        self.failUnlessEqual(s.get(0), s2.get(0))
        self.failUnlessEqual(s.get(1), s2.get(1))
        
