##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2007 Richard Roettger (roettger@in.tum.de)                  #
# Copyright (C) 2008 Dirk Plueger (pflueged@in.tum.de)                      #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
# Copyright (C) 2009 Valeriy Khakhutskyy (khakhutv@in.tum.de)               #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.data.DataEntry import DataEntry
from bin.pysgpp import DataVector
from bin.data.DataContainer import DataContainer

class TestDataContainer(unittest.TestCase):
    
    container = None
    size = None
    dim = None
    vectors = None
    
    #---
    #    It makes following construct: 
    #        [[1,1,1,1,1],  [1, 
    #         [2,2,2,2,2],   2, 
    #         ...            ... 
    #         [42,...,42]]   42]
    #---
    def setUp(self):
         
        self.size = 42
        self.dim = 5
        self.container = DataContainer(self.size,self.dim)
        values = self.container.getValues()
        points = self.container.getPoints()
        self.vectors = []
        for row in xrange(0,self.size):
            vector = DataVector(1,self.dim)
            vector.setAll(row)
            self.vectors.append(vector)
            points.setRow(row,vector)
            values[row] =row
    
    def testNext(self):
        c = 0
        for entry in self.container:
            self.assertEqual(entry.getPoint()[1], self.vectors[c][1])
            self.assertEqual(entry.getValue()[0], c)
            c += 1
        self.assertEqual(c, self.size)

    def testGetTrainDataset(self):
        c = 0
        trainContainer = self.container.getTrainDataset()
        for entry in trainContainer:
            self.assertEqual(entry.getPoint()[1], self.vectors[c][1])
            self.assertEqual(entry.getValue()[0], c)
            c += 1

    def testGetTestDataset(self):
        container = DataContainer(self.container.getPoints(), self.container.getValues(), DataContainer.TEST_CATEGORY)
        c = 0
        testContainer = container.getTestDataset()
        for entry in testContainer:
            self.assertEqual(entry.getPoint()[1], self.vectors[c][1])
            self.assertEqual(entry.getValue()[0], c)
            c += 1

#    def testLoad(self):
#        self.fail("Not implemented")
#
#    def testNormalize(self):
#        self.fail("Not implemented")

    def testCombine(self):
        container = DataContainer(self.container.getPoints(), self.container.getValues(), DataContainer.TEST_CATEGORY)
        self.container = self.container.combine(container)
        self.testGetTrainDataset()
        self.testGetTestDataset()

    def testCreateNullVector(self):
        vector = self.container.createNullVector(self.size, self.dim)
        entry = DataVector(self.dim)
        for row in xrange(self.size):
            vector.getRow(row, entry)
            for index in xrange(self.dim):
                self.assertEqual(entry[index], 0)
        
if __name__=="__main__":
    unittest.main()   