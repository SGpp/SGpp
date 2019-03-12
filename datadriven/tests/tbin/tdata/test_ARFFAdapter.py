# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp.extensions.datadriven.data.ARFFAdapter import ARFFAdapter
from pysgpp import DataVector, DataMatrix


##
# @package tests.tbin.test_ARFFAdapter
# Contains class test_ARFFAdapter::TestARFFAdapter with unittests for @link python.data.ARFFAdapter.ARFFAdapter ARFFAdapter @endlink

##
# Class with unittests for @link python.data.ARFFAdapter.ARFFAdapter ARFFAdapter @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.data.ARFFAdapter.ARFFAdapter ARFFAdapter @endlink
class TestARFFAdapter(unittest.TestCase):


    ##
    # Tests the function @link python.data.ARFFAdapter.ARFFAdapter.save() ARFFAdapter.save() @endlink
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
        points = DataMatrix(size, dim)

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
            self.assertEqual(values[rowIdx], testValues[rowIdx])

        os.remove(filename)


    ##
    # Tests the function @link python.data.ARFFAdapter.ARFFAdapter.loadData() ARFFAdapter.loadData() @endlink
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
        filename = pathlocal + '/../../../datasets/liver/liver-disorders_normalized.arff.gz'
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
            self.assertEqual(values[rowIdx], testValues[rowIdx])


    ##
    # Tests the function @link python.data.ARFFAdapter.ARFFAdapter.loadSpecification() ARFFAdapter.loadSpecification() @endlink
    def testLoadSpecification(self):
        attributes = {
                      "x0":"NUMERIC",
                      "x1":"NUMERIC",
                      "x2":"NUMERIC",
                      "class":"NUMERIC",
                      }
        filename = pathlocal + '/../../../datasets/liver/liver-disorders_normalized.arff.gz'
        adapter = ARFFAdapter(filename)
        spec = adapter.loadSpecification()
        testAttributes = spec.getAttributes()
        self.assertEqual(len(testAttributes), len(attributes))
        for key in testAttributes.keys():
            self.assertEqual(testAttributes[key],attributes[key])

if __name__=="__main__":
    unittest.main()
