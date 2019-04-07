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

from pysgpp.extensions.datadriven.learner.folding.SequentialFoldingPolicy import SequentialFoldingPolicy
from pysgpp.extensions.datadriven.data.DataContainer import DataContainer
from pysgpp import DataVector, DataMatrix


##
# @package tests.tbin.test_SequentialFoldingPolicy
# Contains class test_SequentialFoldingPolicy::TestSequentialFoldingPolicy with unittests for @link python.learner.folding.SequentialFoldingPolicy.SequentialFoldingPolicy SequentialFoldingPolicy @endlink

##
# Class with unittests for @link python.learner.folding.SequentialFoldingPolicy.SequentialFoldingPolicy SequentialFoldingPolicy @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.learner.folding.SequentialFoldingPolicy.SequentialFoldingPolicy SequentialFoldingPolicy @endlink
class TestSequentialFoldingPolicy(unittest.TestCase):


    ## Set up the variables
    def setUp(self):
        self.size = 11
        self.level = 10
        points = DataMatrix(self.size, 1)
        values = DataVector(self.size)
        for i in range(self.size):
            points.set(i, 0, i)
            values[i] = i
        self.dataContainer = DataContainer(points=points, values=values)
        self.policy = SequentialFoldingPolicy(self.dataContainer, self.level)


    ##
    # Tests the function @link python.learner.folding.FoldingPolicy.FoldingPolicy.next() SequentialFoldingPolicy.next() @endlink
    def testNext(self):
        self.assertEqual(self.level, len(self.policy.dataFold))
        for l in self.policy:
            sizeTrain = l.getTrainDataset().getPoints().getSize()
            sizeValidation = l.getTestDataset().getPoints().getSize()
            self.assertEqual(self.size, sizeTrain + sizeValidation)


if __name__=="__main__":
    unittest.main()
