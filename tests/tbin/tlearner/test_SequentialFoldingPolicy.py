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

from bin.learner.folding.SequentialFoldingPolicy import SequentialFoldingPolicy
from bin.data.DataContainer import DataContainer
from pysgpp import DataVector, DataMatrix


##
# @package tests.tbin.test_SequentialFoldingPolicy
# Contains class test_SequentialFoldingPolicy::TestSequentialFoldingPolicy with unittests for @link bin.learner.folding.SequentialFoldingPolicy.SequentialFoldingPolicy SequentialFoldingPolicy @endlink

##
# Class with unittests for @link bin.learner.folding.SequentialFoldingPolicy.SequentialFoldingPolicy SequentialFoldingPolicy @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.learner.folding.SequentialFoldingPolicy.SequentialFoldingPolicy SequentialFoldingPolicy @endlink
class TestSequentialFoldingPolicy(unittest.TestCase):
    
    
    ## Set up the variables
    def setUp(self):
        self.size = 11
        self.level = 10
        points = DataMatrix(self.size, 1)
        values = DataVector(self.size)
        for i in xrange(self.size):
            points.set(i, 0, i)
            values[i] = i
        self.dataContainer = DataContainer(points=points, values=values)
        self.policy = SequentialFoldingPolicy(self.dataContainer, self.level)
    
    
    ##
    # Tests the function @link bin.learner.folding.FoldingPolicy.FoldingPolicy.next() SequentialFoldingPolicy.next() @endlink
    def testNext(self):
        self.assertEqual(self.level, len(self.policy.dataFold))
        for l in self.policy:
            sizeTrain = l.getTrainDataset().getPoints().getSize()
            sizeValidation = l.getTestDataset().getPoints().getSize()
            self.assertEqual(self.size, sizeTrain + sizeValidation)
        
        
if __name__=="__main__":
    unittest.main() 