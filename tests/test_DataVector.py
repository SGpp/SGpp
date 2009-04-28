#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU Lesser General Public License as published  #
# by the Free Software Foundation; either version 3 of the License, or      #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public License  #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################


import unittest

class TestDataVector(unittest.TestCase):
    def testConstructor(self):
        from pysgpp import DataVector
        
        d = DataVector(2)
        self.assertEqual(d.getSize(), 2)
        
        d = DataVector(2,2)
        self.assertEqual(d.getSize(), 2)
        self.assertEqual(d.getDim(), 2)
        self.assertEqual(len(d), 4)

    def testMinMax(self):
        from pysgpp import DataVector
        
        d = DataVector(3)
        for i in xrange(len(d)):
            d[i] = i
            
        self.assertEqual(d.min(0), 0)
        self.assertEqual(d.max(0), 2)
        
        mi, ma = d.minmax(0)
        self.assertEqual(mi, 0)
        self.assertEqual(ma, 2)        

    def testSum(self):
        from pysgpp import DataVector
        
        d = DataVector(3)
        for i in xrange(len(d)):
            d[i] = i
        
        self.assertEqual(d.sum(), 3)

    def testDotProduct(self):
        from pysgpp import DataVector
        
        x = 0
        
        d = DataVector(3)
        for i in xrange(len(d)):
            d[i] = i + 1
            x += d[i] * d[i]
            
        self.assertEqual(d.dotProduct(d), x)

