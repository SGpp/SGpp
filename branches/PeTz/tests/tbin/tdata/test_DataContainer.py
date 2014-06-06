###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de), Dirk Plueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Valeriy Khakhutskyy (khakhutv@in.tum.de)####################################################################

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.data.DataEntry import DataEntry
from bin.pysgpp import DataVector
from bin.data.DataContainer import DataContainer


##
# @package tests.tbin.test_DataContainer
# Contains class test_DataContainer::TestDataContainer with unittests for @link bin.data.DataContainer.DataContainer DataContainer @endlink

##
# Class with unittests for @link bin.data.DataContainer.DataContainer DataContainer @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.data.DataContainer.DataContainer DataContainer @endlink
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
        self.container = DataContainer(self.size,self.dim)
        values = self.container.getValues()
        points = self.container.getPoints()
        self.vectors = []
        for row in xrange(0,self.size):
            vector = DataVector(self.dim)
            vector.setAll(row)
            self.vectors.append(vector)
            points.setRow(row,vector)
            values[row] =row
    
    
    ##
    # Tests the function @link bin.data.DataContainer.DataContainer.next() DataContainer.next() @endlink
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
    # Tests the function @link bin.data.DataContainer.DataContainer.getTrainDataset() DataContainer.getTrainDataset() @endlink
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
    # Tests the function @link bin.data.DataContainer.DataContainer.getTestDataset() DataContainer.getTestDataset() @endlink
    def testGetTestDataset(self):
        container = DataContainer(self.container.getPoints(), self.container.getValues(), DataContainer.TEST_CATEGORY)
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
    # Tests the function @link bin.data.DataContainer.DataContainer.combine() DataContainer.combine() @endlink
    def testCombine(self):
        container = DataContainer(self.container.getPoints(), self.container.getValues(), DataContainer.TEST_CATEGORY)
        self.container = self.container.combine(container)
        self.testGetTrainDataset()
        self.testGetTestDataset()


    ##
    # Tests the function @link bin.data.DataContainer.DataContainer.createNullVector() DataContainer.createNullVector() @endlink
    def testCreateNullVector(self):
        vector = self.container.createNullVector(self.size, self.dim)
        entry = DataVector(self.dim)
        for row in xrange(self.size):
            vector.getRow(row, entry)
            for index in xrange(self.dim):
                self.assertEqual(entry[index], 0)
        
        
        
if __name__=="__main__":
    unittest.main()   