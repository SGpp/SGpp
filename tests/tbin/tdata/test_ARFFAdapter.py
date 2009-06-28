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
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.data.ARFFAdapter import ARFFAdapter
from bin.data.DataVector import DataVector


class TestARFFAdapter(unittest.TestCase):
    

    
    def testSave(self):
        filename = pathlocal + '/datasets/saving.arff.gz'
        testPoints = [[0.307143,0.130137,0.050000],
                      [0.365584,0.105479,0.050000],
                      [0.178571,0.201027,0.050000],
                      [0.272078,0.145548,0.050000],
                      [0.318831,0.065411,0.050000],
                      [0.190260,0.086986,0.050000],
                      [0.190260,0.062329,0.072500],
                      [0.120130,0.068493,0.072500],
                      [0.225325,0.056164,0.072500],
                      [0.213636,0.050000,0.072500]
                     ]
        testValues = [-1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, -1.000000, -1.000000, -1.000000, -1.000000]
        attributes = {
                      "x0":"NUMERIC",
                      "x1":"NUMERIC",
                      "x2":"NUMERIC",
                      "class":"NUMERIC",
                      }
        size = len(testPoints)
        dim = len(testPoints[0])
        point = DataVector(dim)
        points = DataVector(size, dim)
        
        for row in xrange(size):
            for col in xrange(dim):
                point[col] = testPoints[row][col]
            points.setRow(row, point)
            
        adapter = ARFFAdapter(filename)
        adapter.save(points, testValues, attributes)
        
        (points, values) = adapter.loadData().getPointsValues()
        size = len(testPoints)
        dim = len(testPoints[0])
        testVector = DataVector(dim)
        for rowIdx in xrange(size):
            points.getRow(rowIdx, testVector)
            for colIdx in xrange(dim):
                self.assertEqual(testVector[colIdx], testPoints[rowIdx][colIdx])
            self.assertEqual(values[colIdx], testValues[colIdx])
        
        
    def testLoadData(self):
        testPoints = [[0.307143,0.130137,0.050000],
                      [0.365584,0.105479,0.050000],
                      [0.178571,0.201027,0.050000],
                      [0.272078,0.145548,0.050000],
                      [0.318831,0.065411,0.050000],
                      [0.190260,0.086986,0.050000],
                      [0.190260,0.062329,0.072500],
                      [0.120130,0.068493,0.072500],
                      [0.225325,0.056164,0.072500],
                      [0.213636,0.050000,0.072500]
                     ]
        testValues = [-1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, -1.000000, -1.000000, -1.000000, -1.000000]
        filename = pathlocal + '/datasets/liver-disorders_normalized.arff.gz'
        adapter = ARFFAdapter(filename)
        container = adapter.loadData()
        points = container.getPoints()
        values = container.getValues() 
        size = len(testPoints)
        dim = len(testPoints[0])
        testVector = DataVector(dim)
        for rowIdx in xrange(size):
            points.getRow(rowIdx, testVector)
            for colIdx in xrange(dim):
                self.assertEqual(testVector[colIdx], testPoints[rowIdx][colIdx])
            self.assertEqual(values[colIdx], testValues[colIdx])
        

    def testLoadSpecification(self):
        attributes = {
                      "x0":"NUMERIC",
                      "x1":"NUMERIC",
                      "x2":"NUMERIC",
                      "class":"NUMERIC",
                      }
        filename = pathlocal + '/datasets/liver-disorders_normalized.arff.gz'
        adapter = ARFFAdapter(filename)        
        spec = adapter.loadSpecification()
        self.assertEqual(len(spec), len(attributes))
        for key in spec.keys():
            self.assertEqual(spec[key],attributes[key])
        
if __name__=="__main__":
    unittest.main() 