# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp.extensions.datadriven.data.DataEntry import DataEntry
from pysgpp import DataVector
from pysgpp.extensions.datadriven.data.DataContainer import DataContainer


##
# @package tests.tbin.test_DataContainer
# Contains class test_DataContainer::TestDataContainer with unittests for @link python.data.DataContainer.DataContainer DataContainer @endlink

##
# Class with unittests for @link python.data.DataContainer.DataContainer DataContainer @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.data.DataContainer.DataContainer DataContainer @endlink
class TestDataContainer(unittest.TestCase):


    ##
    # Set up the variables
    #    It makes following construct:
    #        [[1,1,1,1,1],  [1,
    #         [2,2,2,2,2],   2,
    #         ...            ...
    #         [42,...,42]]   42]
    def setUp(self):
        self.size = 42
        self.dim = 5
        self.container = DataContainer(size=self.size,dim=self.dim)
        values = self.container.getValues()
        points = self.container.getPoints()
        self.vectors = []
        for row in range(0,self.size):
            vector = DataVector(self.dim)
            vector.setAll(row)
            self.vectors.append(vector)
            points.setRow(row,vector)
            values[row] =row


    ##
    # Tests the function @link python.data.DataContainer.DataContainer.next() DataContainer.next() @endlink
    def testNext(self):
        c = 0
        for entry in self.container:
            point = entry.getPoint()
            point.sub(self.vectors[c])
            self.assertEqual(point.l2Norm(), 0)
            self.assertEqual(entry.getValue(), c)
            c += 1
        self.assertEqual(c, self.size)


    ##
    # Tests the function @link python.data.DataContainer.DataContainer.getTrainDataset() DataContainer.getTrainDataset() @endlink
    def testGetTrainDataset(self):
        c = 0
        trainContainer = self.container.getTrainDataset()
        for entry in trainContainer:
            point = entry.getPoint()
            point.sub(self.vectors[c])
            self.assertEqual(point.l2Norm(), 0)
            self.assertEqual(entry.getValue(), c)
            c += 1


    ##
    # Tests the function @link python.data.DataContainer.DataContainer.getTestDataset() DataContainer.getTestDataset() @endlink
    def testGetTestDataset(self):
        container = DataContainer(points=self.container.getPoints(), values=self.container.getValues(), name=DataContainer.TEST_CATEGORY)
        c = 0
        testContainer = container.getTestDataset()
        for entry in testContainer:
            point = entry.getPoint()
            point.sub(self.vectors[c])
            self.assertEqual(point.l2Norm(), 0)
            self.assertEqual(entry.getValue(), c)
            c += 1

#    def testLoad(self):
#        self.fail("Not implemented")
#
#    def testNormalize(self):
#        self.fail("Not implemented")


    ##
    # Tests the function @link python.data.DataContainer.DataContainer.combine() DataContainer.combine() @endlink
    def testCombine(self):
        container = DataContainer(points=self.container.getPoints(), values=self.container.getValues(), name=DataContainer.TEST_CATEGORY)
        self.container = self.container.combine(container)
        self.testGetTrainDataset()
        self.testGetTestDataset()


    ##
    # Tests the function @link python.data.DataContainer.DataContainer.createNullVector() DataContainer.createNullVector() @endlink
    def testCreateNullVector(self):
        vector = self.container.createNullVector(self.size, self.dim)
        entry = DataVector(self.dim)
        for row in range(self.size):
            vector.getRow(row, entry)
            for index in range(self.dim):
                self.assertEqual(entry[index], 0)



if __name__=="__main__":
    unittest.main()
