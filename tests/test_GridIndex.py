###############################################################################
# Copyright (C) 2014 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
# @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

import math
from pysgpp import GridIndex
import unittest

def ccKnot(l, i):
    """Return Clenshaw-Curtis knot with given level and index."""
    return 0.5 * (math.cos(math.pi * (1.0 - i / 2.0**l)) + 1.0)

class TestSGridIndex(unittest.TestCase):
    def testConstructor(self):
        """Tests copy constructor"""
        s = GridIndex(2)
        s.set(0,1,1)
        s.set(0,2,3)
        
        s2 = GridIndex(s)
        self.assertEqual(s.get(0), s2.get(0))
    
    def testAssign(self):
        """Tests assignment"""
        s = GridIndex(2)
        s.set(0,1,1)
        s.set(1,2,3)
        
        s2 = GridIndex(5)
        
        s2.assign(s)
        self.assertEqual(s.get(0), s2.get(0))
        self.assertEqual(s.get(1), s2.get(1))
    
    def testClenshawCurtis(self):
        """Tests Clenshaw-Curtis coordinates"""
        d = 5
        s = GridIndex(d)
        
        # set grid point levels/indices
        s.set(0, 0, 0)
        s.set(1, 0, 1)
        s.set(2, 1, 1)
        s.set(3, 3, 5)
        s.set(4, 4, 11)
        s.setPointDistribution(GridIndex.ClenshawCurtis)
        
        for t in range(d):
            # check twice to test coordinate caching in sg::base::GridIndex
            self.assertAlmostEqual(ccKnot(s.getLevel(t), s.getIndex(t)), s.abs(t))
            self.assertAlmostEqual(ccKnot(s.getLevel(t), s.getIndex(t)), s.abs(t))
