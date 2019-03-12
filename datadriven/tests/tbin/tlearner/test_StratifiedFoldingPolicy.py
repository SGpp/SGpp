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

from pysgpp.extensions.datadriven.learner.folding.StratifiedFoldingPolicy import StratifiedFoldingPolicy
from pysgpp.extensions.datadriven.data.DataContainer import DataContainer
from pysgpp import DataVector, DataMatrix


##
# @package tests.tbin.test_StratifiedFoldingPolicy
# Contains class test_StratifiedFoldingPolicy::TestStratifiedFoldingPolicy with unittests for @link python.learner.folding.StratifiedFoldingPolicy.StratifiedFoldingPolicy StratifiedFoldingPolicy @endlink

##
# Class with unittests for @link python.learner.folding.StratifiedFoldingPolicy.StratifiedFoldingPolicy StratifiedFoldingPolicy @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.learner.folding.StratifiedFoldingPolicy.StratifiedFoldingPolicy StratifiedFoldingPolicy @endlink
class TestStratifiedFoldingPolicy(unittest.TestCase):


    ## Set up the variables
    def setUp(self):
        self.size = 9
        self.level = 4
        points = DataMatrix(self.size, 1)
        values = DataVector(self.size)
        for i in range(self.size):
            points.set(i, 0, i)
            values[i] = -1 if i < self.size//2 else 1
        self.dataContainer = DataContainer(points=points, values=values)
        self.policy = StratifiedFoldingPolicy(self.dataContainer, self.level)


    ##
    # Tests the function @link python.learner.folding.FoldingPolicy.FoldingPolicy.next() StratifiedFoldingPolicy.next() @endlink
    def testNext(self):
        validationCorrectData = [[4,0],[5,1], [6,2], [7,8,3]]
        self.assertEqual(self.level, len(self.policy.dataFold))
        i = 0
        for l in self.policy:
            trainPoints = l.getTrainDataset().getPoints()
            validationPoints = l.getTestDataset().getPoints()
            validCorrectD = set(validationCorrectData[i])
            trainCorrectD = set(range(self.size)) - validCorrectD

            self.assertEqual(trainPoints.getNrows(), len(trainCorrectD))
            self.assertEqual(validationPoints.getSize(), len(validCorrectD))
            for k in range(trainPoints.getNrows()):
                self.assertTrue(trainPoints.get(k,0)  in trainCorrectD)
            for k in range(validationPoints.getSize()):
                self.assertTrue(validationPoints.get(k,0) in validCorrectD)

            i += 1




if __name__=="__main__":
    unittest.main()
